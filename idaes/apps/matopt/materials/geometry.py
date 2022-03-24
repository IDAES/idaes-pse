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
__all__ = [
    "Parallelepiped",
    "RectPrism",
    "Cube",
    "Rhombohedron",
    "Cuboctahedron",
    "Cylinder",
    "CylindricalSector",
]

from abc import abstractmethod, ABC
from copy import deepcopy
from math import cos, sin, sqrt

import numpy as np

from .transform_func import ShiftFunc, ScaleFunc, RotateFunc, ReflectFunc
from ..util.util import areEqual


class Shape(object):
    """ """

    DBL_TOL = 1e-5
    DEFAULT_ALIGNMENT = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)

    # === STANDARD CONSTRUCTOR
    def __init__(self, Anchor, Alignment=None):
        self._Anchor = Anchor
        self._Alignment = Shape.DEFAULT_ALIGNMENT if Alignment is None else Alignment

    # === ASSERTION OF CLASS DESIGN
    def isConsistentWithDesign(self):
        """ """
        if self.Anchor is None:
            print("alpha")
            return False
        if type(self.Alignment) is not np.ndarray:
            print("A")
            return False
        if self.Alignment.shape != (3, 3):
            print("B")
            return False
        if not areEqual(np.linalg.det(self.Alignment), 1.0, Shape.DBL_TOL):
            print("C")
            return False
        for AlignmentAxis in self.Alignment:
            if not areEqual(np.linalg.norm(AlignmentAxis), 1.0, Shape.DBL_TOL):
                print("D")
                return False
        return True

    # === MANIPULATION METHODS
    def applyTransF(self, TransF):
        """

        Args:
            TransF:

        Returns:

        """
        if type(TransF) is ShiftFunc or type(TransF) is ScaleFunc:
            TransF.transform(self._Anchor)
        elif type(TransF) is RotateFunc or type(TransF) is ReflectFunc:
            TransF.transform(self._Anchor)
            TransF.transform(self._Alignment)
            assert self.isConsistentWithDesign()
        else:
            raise TypeError
        assert self.isConsistentWithDesign()

    def shift(self, Shift):
        """

        Args:
            Shift:

        Returns:

        """
        if type(Shift) is ShiftFunc:
            self.applyTransF(Shift)
        elif type(Shift) is np.ndarray:
            self.applyTransF(ShiftFunc(Shift))
        else:
            raise TypeError

    def scale(self, Scale, OriginOfScale=None):
        """

        Args:
            Scale: param OriginOfScale:  (Default value = None)
            OriginOfScale:  (Default value = None)

        Returns:

        """
        if type(Scale) is ScaleFunc:
            self.applyTransF(Scale)
        elif type(Scale) is np.ndarray:
            self.applyTransF(ScaleFunc(Scale, OriginOfScale))
        elif type(Scale) is float or type(Scale) is int:
            self.applyTransF(ScaleFunc(np.array([Scale, Scale, Scale], dtype=float)))
        else:
            raise TypeError

    def rotate(self, Rotation, OriginOfRotation=None):
        """

        Args:
            Rotation: param OriginOfRotation:  (Default value = None)
            OriginOfRotation:  (Default value = None)

        Returns:

        """
        if type(Rotation) is RotateFunc:
            self.applyTransF(Rotation)
        elif type(Rotation) is np.ndarray:
            self.applyTransF(RotateFunc(Rotation, OriginOfRotation))
        else:
            raise TypeError

    def reflect(self, Reflection):
        """

        Args:
            Reflection:

        Returns:

        """
        if type(Reflection) is ReflectFunc:
            self.applyTransF(Reflection)
        else:
            raise TypeError

    # === PROPERTY EVALUATION METHODS
    @abstractmethod
    def isInShape(self, P):
        """

        Args:
            P:

        Returns:

        """
        raise NotImplementedError

    __contains__ = isInShape

    @abstractmethod
    def getBounds(self):
        """ """
        raise NotImplementedError

    # === BASIC QUERY METHODS
    @property
    def Anchor(self):
        """ """
        return self._Anchor

    @property
    def Alignment(self):
        """ """
        return self._Alignment


class Polyhedron(Shape, ABC):
    """ """

    # === STANDARD CONSTRUCTOR
    def __init__(self, V, F, Anchor, Alignment=None):
        Shape.__init__(self, Anchor, Alignment)
        self._V = V
        self._F = F
        self._FacetNorms = self.__calcFacetNorms()
        self._FacetDirections = self.__calcFacetDirections()
        assert self.isConsistentWithDesign()

    # === ASSERTION OF CLASS DESIGN
    def isConsistentWithDesign(self):
        """ """
        if type(self.V) is not list:
            return False
        if len(self.V) < 4:
            return False
        if type(self.F) is not list:
            return False
        if len(self.F) < 4:
            return False
        for Facet in self.F:
            if type(Facet) is not list:
                return False
            if len(Facet) < 3:
                return False
        return Shape.isConsistentWithDesign(self)

    # === AUXILIARY METHODS
    def __calcFacetNorms(self):
        FacetNorms = []
        for Facet in self.F:
            Norm = np.cross(
                self.V[Facet[1]] - self.V[Facet[0]], self.V[Facet[2]] - self.V[Facet[0]]
            )
            Norm /= np.linalg.norm(Norm)
            FacetNorms.append(Norm)
        return FacetNorms

    def __calcFacetDirections(self):
        FacetDirections = []
        for Facet in self.F:
            assert len(Facet) >= 3
            FacetDirection = np.array([0, 0, 0], dtype=float)
            for v in Facet:
                FacetDirection += self.V[v]
            FacetDirection /= len(Facet)
            FacetDirection -= self.Anchor
            FacetDirections.append(FacetDirection)
        return FacetDirections

    # === MANIPULATION METHODS
    def applyTransF(self, TransF):
        """

        Args:
            TransF:

        Returns:

        """
        if type(TransF) is ShiftFunc:
            for v in self._V:
                TransF.transform(v)
        elif type(TransF) is ScaleFunc:
            for v in self._V:
                TransF.transform(v)
            for direction in self._FacetDirections:
                TransF.transform(direction)
            for norm in self._FacetNorms:
                norm /= TransF.Scale
        elif type(TransF) is RotateFunc:
            for v in self._V:
                TransF.transform(v)
            for direction in self._FacetDirections:
                TransF.transformDirection(direction)
            for norm in self._FacetNorms:
                TransF.transformDirection(norm)
        else:
            raise TypeError(
                "MatOpt does not support this transformation type. Please contact MatOpt developer for "
                "potential feature addition."
            )
        Shape.applyTransF(self, TransF)
        assert self.isConsistentWithDesign()

    # === PROPERTY EVALUATION METHODS
    def isInShape(self, P, tol=Shape.DBL_TOL):
        """

        Args:
            P: param tol:  (Default value = Shape.DBL_TOL)
            tol:  (Default value = Shape.DBL_TOL)

        Returns:

        """
        for f in range(len(self.F)):
            if not self.satisfiesFacet(P, f, tol):
                return False
        return True

    __contains__ = isInShape

    def satisfiesFacet(self, P, f, tol=Shape.DBL_TOL):
        """

        Args:
            P: param f:
            tol: Default value = Shape.DBL_TOL)
            f:

        Returns:

        """
        P01 = P - self.V[self.F[f][0]]
        return np.inner(P01, self.FacetNorms[f]) < tol

    def getBounds(self):
        """ """
        MinP = deepcopy(self.V[0])
        MaxP = deepcopy(self.V[0])
        for v in self.V:
            MinP = np.minimum(MinP, v)
            MaxP = np.maximum(MaxP, v)
        return MinP, MaxP

    # === BASIC QUERY METHODS
    @property
    def V(self):
        """ """
        return self._V

    @property
    def F(self):
        """ """
        return self._F

    @property
    def FacetNorms(self):
        """ """
        return self._FacetNorms

    @property
    def FacetDirections(self):
        """ """
        return self._FacetDirections


class Cuboctahedron(Polyhedron, ABC):
    """ """

    # === STANDARD CONSTRUCTOR
    def __init__(self, R, Center=None):
        self._R = sqrt(2)
        V = [
            np.array([1, 0, 1], dtype=float),
            np.array([0, -1, 1], dtype=float),
            np.array([-1, 0, 1], dtype=float),
            np.array([0, 1, 1], dtype=float),
            np.array([1, 1, 0], dtype=float),
            np.array([1, -1, 0], dtype=float),
            np.array([-1, -1, 0], dtype=float),
            np.array([-1, 1, 0], dtype=float),
            np.array([1, 0, -1], dtype=float),
            np.array([0, -1, -1], dtype=float),
            np.array([-1, 0, -1], dtype=float),
            np.array([0, 1, -1], dtype=float),
        ]
        F = [
            [0, 3, 2, 1],
            [1, 6, 9, 5],
            [2, 7, 10, 6],
            [3, 4, 11, 7],
            [0, 5, 8, 4],
            [8, 9, 10, 11],
            [1, 2, 6],
            [2, 3, 7],
            [0, 4, 3],
            [0, 1, 5],
            [6, 10, 9],
            [7, 11, 10],
            [4, 8, 11],
            [5, 9, 8],
        ]
        Polyhedron.__init__(self, V, F, np.array([0, 0, 0], dtype=float))
        self.scale(R / sqrt(2))
        if Center is not None:
            self.shift(Center)

    # === MANIPULATION METHODS
    def applyTransF(self, TransF):
        """

        Args:
            TransF:

        Returns:

        """
        if isinstance(TransF, ScaleFunc):
            if TransF.isIsometric:
                self._R *= TransF.Scale[0]
            else:
                raise ValueError(
                    "Cuboctahedron applyTransF: Can only scale isometrically"
                )
        Polyhedron.applyTransF(self, TransF)
        assert self.isConsistentWithDesign()


class Parallelepiped(Polyhedron, ABC):
    """ """

    # === STANDARD CONSTRUCTOR
    def __init__(self, Vx, Vy, Vz, BotBackLeftCorner=None):
        self._Vx = Vx
        self._Vy = Vy
        self._Vz = Vz
        if np.inner(np.cross(Vx, Vy), Vz) < Parallelepiped.DBL_TOL:
            raise ValueError("Vx,Vy,Vz must have a positive box product.")
        if BotBackLeftCorner is None:
            BotBackLeftCorner = np.array([0, 0, 0], dtype=float)

        V = [
            BotBackLeftCorner,
            BotBackLeftCorner + Vx,
            BotBackLeftCorner + Vy,
            BotBackLeftCorner + Vz,
            BotBackLeftCorner + Vx + Vy,
            BotBackLeftCorner + Vx + Vz,
            BotBackLeftCorner + Vy + Vz,
            BotBackLeftCorner + Vx + Vy + Vz,
        ]
        F = [
            [0, 1, 5, 3],
            [0, 3, 6, 2],
            [2, 6, 7, 4],
            [1, 4, 7, 5],
            [3, 5, 7, 6],
            [0, 2, 4, 1],
        ]
        Polyhedron.__init__(self, V, F, Anchor=deepcopy(BotBackLeftCorner))

    # === CONSTRUCTOR - From edge lengths and angles
    @classmethod
    def fromEdgesAndAngles(cls, A, B, C, alpha, beta, gamma, BotBackLeftCorner=None):
        """

        Args:
            A: param B:
            C: param alpha:
            beta: param gamma:
            BotBackLeftCorner: Default value = None)
            B:
            alpha:
            gamma:

        Returns:

        """
        Vx = np.array([A, 0, 0])
        Vy = np.array([B * np.cos(gamma), B * np.sin(gamma), 0])
        zcomp = np.sqrt(1 - pow(np.cos(alpha), 2) - pow(np.cos(beta), 2))
        Vz = np.array([C * np.cos(alpha), C * np.cos(beta), C * zcomp])
        return cls(Vx, Vy, Vz, BotBackLeftCorner)

    # === CONSTRUCTOR - From POSCAR file header information
    @classmethod
    def fromPOSCAR(cls, filename):
        """

        Args:
            filename:

        Returns:

        """
        with open(filename, "r") as infile:
            CommentLine = infile.readline()
            GSLine = infile.readline().split()
            GS = float(GSLine[0])
            VxLine = infile.readline().split()
            Vx = GS * np.array(
                [float(VxLine[0]), float(VxLine[1]), float(VxLine[2])], dtype=float
            )
            VyLine = infile.readline().split()
            Vy = GS * np.array(
                [float(VyLine[0]), float(VyLine[1]), float(VyLine[2])], dtype=float
            )
            VzLine = infile.readline().split()
            Vz = GS * np.array(
                [float(VzLine[0]), float(VzLine[1]), float(VzLine[2])], dtype=float
            )
        return cls(Vx, Vy, Vz)

    # === MANIPULATION METHODS
    def applyTransF(self, TransF):
        """

        Args:
            TransF:

        Returns:

        """
        if isinstance(TransF, ScaleFunc):
            self._Vx *= TransF.Scale[0]
            self._Vy *= TransF.Scale[1]
            self._Vz *= TransF.Scale[2]
        Polyhedron.applyTransF(self, TransF)
        assert self.isConsistentWithDesign()

        # === PROPERTY EVALUATION METHODS

    def getFractionalCoords(self, P, blnRelativeToCenter=False):
        """

        Args:
            P: param blnRelativeToCenter:  (Default value = False)
            blnRelativeToCenter:  (Default value = False)

        Returns:

        """
        FracOrigin = self.getCenter() if blnRelativeToCenter else self.Anchor
        Omega = self.getVolume()
        fracX = 1.0 / Omega * np.inner(P - FracOrigin, np.cross(self.Vy, self.Vz))
        fracY = 1.0 / Omega * np.inner(P - FracOrigin, np.cross(self.Vz, self.Vx))
        fracZ = 1.0 / Omega * np.inner(P - FracOrigin, np.cross(self.Vx, self.Vy))
        return np.array([fracX, fracY, fracZ])

    def getCenter(self):
        """ """
        return sum(self.V) / len(self.V)

    def getVolume(self):
        """ """
        return np.dot(np.cross(self.Vx, self.Vy), self.Vz)

    # === BASIC QUERY METHODS
    @property
    def Vx(self):
        """ """
        return self._Vx

    @property
    def Vy(self):
        """ """
        return self._Vy

    @property
    def Vz(self):
        """ """
        return self._Vz

    @property
    def isUnit(self):
        """ """
        if not areEqual(np.linalg.norm(self.Vx), 1.0, Shape.DBL_TOL):
            return False
        if not areEqual(np.linalg.norm(self.Vy), 1.0, Shape.DBL_TOL):
            return False
        if not areEqual(np.linalg.norm(self.Vz), 1.0, Shape.DBL_TOL):
            return False
        return True


class Rhombohedron(Parallelepiped, ABC):
    """ """

    # === STANDARD CONSTRUCTOR
    def __init__(self, L, Alpha, BotBackLeftCorner=None):
        self._L = L
        self._Alpha = Alpha
        Lx = np.array([L, 0, 0])
        Ly = np.array([L * cos(Alpha), L * sin(Alpha), 0])
        Lz = np.array(
            [
                L * cos(Alpha),
                L * (cos(Alpha) - cos(Alpha) ** 2) / sin(Alpha),
                L * sqrt(1 - 3 * cos(Alpha) ** 2 + 2 * cos(Alpha) ** 3) / sin(Alpha),
            ]
        )
        Parallelepiped.__init__(self, Lx, Ly, Lz, BotBackLeftCorner)

    # === MANIPULATION METHODS
    def applyTransF(self, TransF):
        """

        Args:
            TransF:

        Returns:

        """
        if isinstance(TransF, ScaleFunc):
            if TransF.isIsometric:
                self._L *= TransF.Scale[0]
            else:
                raise ValueError(
                    "Rhombohedron applyTransF: Can only scale isometrically"
                )
        Polyhedron.applyTransF(self, TransF)
        assert self.isConsistentWithDesign()

        # === PROPERTY EVALUATION METHODS

    @property
    def L(self):
        """ """
        return self._L

    @property
    def Alpha(self):
        """ """
        return self._Alpha


class RectPrism(Parallelepiped, ABC):
    """ """

    # === STANDARD CONSTRUCTOR
    def __init__(self, Lx, Ly, Lz, BotBackLeftCorner=None):
        self._Lx = Lx
        self._Ly = Ly
        self._Lz = Lz
        Parallelepiped.__init__(
            self,
            np.array([Lx, 0, 0], dtype=float),
            np.array([0, Ly, 0], dtype=float),
            np.array([0, 0, Lz], dtype=float),
            BotBackLeftCorner,
        )

    # === CONSTRUCTOR - Bounding box of a list of points
    @classmethod
    def fromPointsBBox(cls, Pts):
        """

        Args:
            Pts:

        Returns:

        """
        # import code; code.interact(local=dict(locals(),**globals()));
        minX = min(_[0] for _ in Pts)
        maxX = max(_[0] for _ in Pts)
        minY = min(_[1] for _ in Pts)
        maxY = max(_[1] for _ in Pts)
        minZ = min(_[2] for _ in Pts)
        maxZ = max(_[2] for _ in Pts)
        return cls(maxX - minX, maxY - minY, maxZ - minZ, np.array([minX, minY, minZ]))

    # === MANIPULATION METHODS
    def applyTransF(self, TransF):
        """

        Args:
            TransF:

        Returns:

        """
        if isinstance(TransF, ScaleFunc):
            self._Lx *= TransF.Scale[0]
            self._Ly *= TransF.Scale[1]
            self._Lz *= TransF.Scale[2]
        Parallelepiped.applyTransF(self, TransF)
        assert self.isConsistentWithDesign()

        # === BASIC QUERY METHODS

    @property
    def Lx(self):
        """ """
        return self._Lx

    @property
    def Ly(self):
        """ """
        return self._Ly

    @property
    def Lz(self):
        """ """
        return self._Lz


class Cube(RectPrism, ABC):
    """ """

    # === STANDARD CONSTRUCTOR
    def __init__(self, L, BotBackLeftCorner=None):
        self._L = L
        RectPrism.__init__(self, L, L, L, BotBackLeftCorner)

    # === MANIPULATION METHODS
    def applyTransF(self, TransF):
        """

        Args:
            TransF:

        Returns:

        """
        if isinstance(TransF, ScaleFunc):
            if TransF.isIsometric:
                self._L *= TransF.Scale[0]
            else:
                raise ValueError("Cube applyTransF: Can only scale isometrically")
        RectPrism.applyTransF(self, TransF)
        assert self.isConsistentWithDesign()

        # === BASIC QUERY METHODS

    @property
    def L(self):
        """ """
        return self._L


class Cylinder(Shape, ABC):
    """Object to create, manipulate and utilize Cylinder shapes."""

    # === STANDARD CONSTRUCTOR
    def __init__(self, Po, R, H, Vh=np.array([0, 0, 1], dtype=float), Alignment=None):
        Shape.__init__(self, Po, Alignment)
        self._R = R
        self._H = H
        self._Vh = Vh * self._H / np.linalg.norm(Vh)
        assert self.isConsistentWithDesign()

    # === ASSERTION OF CLASS DESIGN
    def isConsistentWithDesign(self):
        if len(self.Vh) != 3:
            return False
        return self.R > 0 and self.H > 0 and Shape.isConsistentWithDesign(self)

    # === MANIPULATION METHODS
    def applyTransF(self, TransF):
        """Apply the input transform function.

        Args:
            TransF: Transform function.
        """
        if type(TransF) is ScaleFunc:
            if areEqual(TransF.Scale[0], TransF.Scale[1], Shape.DBL_TOL):
                self._R *= TransF.Scale[0]
                self._H *= TransF.Scale[2]
                TransF.transform(self._Vh)
            else:
                raise ValueError(
                    "The scaling ratio along the radius directions are not consistent"
                )
        elif (type(TransF) is RotateFunc) or (type(TransF) is ReflectFunc):
            TransF.transform(self._Vh)
        elif not (type(TransF) is ShiftFunc):
            raise TypeError(
                "MatOpt does not support this transformation type. Please contact MatOpt developer for "
                "potential feature addition."
            )
        Shape.applyTransF(self, TransF)
        assert self.isConsistentWithDesign()

    # === PROPERTY EVALUATION METHODS
    def isInShape(self, P, tol=Shape.DBL_TOL):
        """Check if a point is inside the Cylinder shape within certain tolerance.

        Args:
            P: Input point
            tol: Tolerance (Default value = Shape.DBL_TOL)

        Returns:
            True/False
        """
        OP = P - self.Anchor
        projH = np.dot(OP, self.Vh) / self.H
        if (projH > -tol) and (projH < self.H + tol):
            normOP = np.linalg.norm(OP)
            projR = sqrt(normOP**2 - projH**2 + tol)
            return projR < self.R + tol
        else:
            return False

    __contains__ = isInShape

    def getBounds(self):
        """Get the bounds of each coordinate in the form of an artificial point."""
        MinP = (
            np.minimum(np.zeros(3, dtype=float), self.Vh)
            + self.Anchor
            - np.ones(3, dtype=float) * self.R
        )
        MaxP = (
            np.maximum(np.zeros(3, dtype=float), self.Vh)
            + self.Anchor
            + np.ones(3, dtype=float) * self.R
        )
        return MinP, MaxP

    # === BASIC QUERY METHODS
    @property
    def R(self):
        """Radius."""
        return self._R

    @property
    def H(self):
        """Height."""
        return self._H

    @property
    def Vh(self):
        """Center axis vector."""
        return self._Vh


class CylindricalSector(Shape, ABC):
    """Object to create, manipulate and utilize Cylinder shapes."""

    # === STANDARD CONSTRUCTOR
    def __init__(self, Po, R, H, Va, Vb, Vh, Alignment=None):
        Shape.__init__(self, Po, Alignment)
        self._R = R
        self._H = H
        self._Vh = Vh * self._H / np.linalg.norm(Vh)
        self._Va = Va * self._R / np.linalg.norm(Va)
        self._Vb = Vb * self._R / np.linalg.norm(Vb)
        self._Norms = self.setNorms()
        assert self.isConsistentWithDesign()

    # === ASSERTION OF CLASS DESIGN
    def isConsistentWithDesign(self):
        if (
            self.R > 0
            and self.H > 0
            and len(self.Vh) == 3
            and len(self.Va) == 3
            and len(self.Vb) == 3
            and len(self.Norms) == 3
        ):
            if areEqual(np.dot(self.Va, self.Vh), 0.0, Shape.DBL_TOL) and areEqual(
                np.dot(self.Vb, self.Vh), 0.0, Shape.DBL_TOL
            ):
                return Shape.isConsistentWithDesign(self)
        return False

    # === AUXILIARY METHODS
    def setNorms(self):
        return [
            np.cross(self._Vb, self._Va),
            np.cross(self._Va, self._Vh),
            np.cross(self._Vh, self._Vb),
        ]

    # === MANIPULATION METHODS
    def applyTransF(self, TransF):
        """Apply the input transform function.

        Args:
            TransF: Transform function.
        """
        if type(TransF) is ScaleFunc:
            if areEqual(TransF.Scale[0], TransF.Scale[1], Shape.DBL_TOL):
                self._R *= TransF.Scale[0]
                self._H *= TransF.Scale[2]
                TransF.transform(self._Vh)
            else:
                raise ValueError(
                    "The scaling ratio along the radius directions are not consistent"
                )
        elif (type(TransF) is RotateFunc) or (type(TransF) is ReflectFunc):
            TransF.transform(self._Va)
            TransF.transform(self._Vb)
            TransF.transform(self._Vh)
            self._Norms = self.setNorms()
        elif not (type(TransF) is ShiftFunc):
            raise TypeError(
                "MatOpt does not support this transformation type. Please contact MatOpt developer for "
                "potential feature addition."
            )
        Shape.applyTransF(self, TransF)
        assert self.isConsistentWithDesign()

    # === PROPERTY EVALUATION METHODS
    def isInShape(self, P, tol=Shape.DBL_TOL):
        """Check if a point is inside the Cylinder shape within certain tolerance.

        Args:
            P: Input point
            tol: Tolerance (Default value = Shape.DBL_TOL)

        Returns:
            True/False
        """
        OP = P - self.Anchor
        for n in self.Norms:
            if np.dot(n, OP) > tol:
                return False
        if np.dot(-self.Norms[0], OP - self.Vh) > tol:
            return False
        projH = np.linalg.norm(np.dot(OP, self.Vh)) / self.H
        normP = np.linalg.norm(OP)
        if sqrt(normP**2 - projH**2 + tol) > self.R + tol:
            return False
        return True

    __contains__ = isInShape

    # === BASIC QUERY METHODS
    @property
    def R(self):
        """Radius."""
        return self._R

    @property
    def H(self):
        """Height."""
        return self._H

    @property
    def Vh(self):
        """Center axis vector."""
        return self._Vh

    @property
    def Va(self):
        """Left edge vector."""
        return self._Va

    @property
    def Vb(self):
        """Right edge vector."""
        return self._Vb

    @property
    def Norms(self):
        """Norms."""
        return self._Norms
