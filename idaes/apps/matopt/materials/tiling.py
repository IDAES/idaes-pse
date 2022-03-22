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
from abc import abstractmethod, ABC
import numpy as np

from .geometry import Parallelepiped, Cylinder, CylindricalSector
from .design import Design


class Tiling(object):
    """ """

    def __init__(self):
        pass

    # === PROPERTY EVALUATION METHODS
    @abstractmethod
    def transformInsideTile(self, P):
        """Transform point to lie inside a tile.

        Args:
            P (numpy.ndarray): Point to modify to be inside tile.

        Returns:
            None.

        """
        raise NotImplementedError

    @abstractmethod
    def replicateDesign(self, D, nTiles, OldToNewIndices=None, AuxPropMap=None):
        """Create a larger Design by tiling a smaller one.

        Args:
            D (Design): A smaller Design to tile.
            nTiles (int/numpy.ndarray): A specifier for number of
                tiles to replicate in periodic directions.
            OldToNewIndices (dict<int,list<int>>): Optional, a
                dictionary to store the mapping of tiled points
                indices from the smaller Design to the larger.
                (Default value = None)
            AuxPropMap (dict<tuple<string,string>,dict<int,float>>):
                Optional, a mapping of locations to properties.
                NOTE: It is modified in place, so make a copy before.
                (Default value = None)

        Returns:
            (Design): A larger, periodic Design

        """
        raise NotImplementedError


class LinearTiling(Tiling, ABC):
    """Class to manage linear periodicity."""

    DBL_TOL = 1e-5

    # === STANDARD CONSTRUCTOR
    def __init__(self, TilingDirections_, argShape=None):
        super().__init__()
        self._TileShape = argShape
        self._TilingDirections = TilingDirections_

    # === CONSTRUCTOR - From Parallelepiped
    @classmethod
    def fromParallelepiped(cls, argShape):
        assert (
            type(argShape) == Parallelepiped
        ), "The input shape is not an instance of Parallelepiped."
        TilingDirections_ = [argShape.Vz, -argShape.Vz]
        return cls(TilingDirections_, argShape)

    # === CONSTRUCTOR - From Cylinder or CylindricalSector
    @classmethod
    def fromCylindricalShape(cls, argShape):
        assert (
            type(argShape) == Cylinder or CylindricalSector
        ), "The input shape is not an instance of Cylinder or CylindricalSector."
        TilingDirections_ = [argShape.Vh, -argShape.Vh]
        return cls(TilingDirections_, argShape)

    # === CONSTRUCTOR - From POSCAR files
    @classmethod
    def fromPOSCAR(cls, filename):
        return cls.fromParallelepiped(Parallelepiped.fromPOSCAR(filename))

    # === BASIC QUERY METHODS
    @property
    def TileShape(self):
        return self._TileShape

    @property
    def TilingDirections(self):
        return self._TilingDirections

    @property
    def V(self):
        return self._TilingDirections[0]


class PlanarTiling(Tiling):
    """ """

    DBL_TOL = 1e-5

    # === STANDARD CONSTRUCTOR
    def __init__(self, Parallelepiped_):
        super().__init__()
        self._TileShape = Parallelepiped_
        self._TilingDirections = []
        self._TilingDirections.append(-Parallelepiped_.Vx)
        self._TilingDirections.append(-Parallelepiped_.Vy)
        self._TilingDirections.append(Parallelepiped_.Vx)
        self._TilingDirections.append(Parallelepiped_.Vy)
        self._TilingDirections.append(-Parallelepiped_.Vx - Parallelepiped_.Vy)
        self._TilingDirections.append(-Parallelepiped_.Vx + Parallelepiped_.Vy)
        self._TilingDirections.append(Parallelepiped_.Vx - Parallelepiped_.Vy)
        self._TilingDirections.append(Parallelepiped_.Vx + Parallelepiped_.Vy)

    # === CONSTRUCTOR - From POSCAR file
    @classmethod
    def fromPOSCAR(cls, filename):
        """

        Args:
            filename:

        Returns:

        """
        return cls(Parallelepiped.fromPOSCAR(filename))

    # === PROPERTY EVALUATION METHODS
    def transformInsideTile(self, P, EdgeTol=DBL_TOL):
        """

        Args:
            P: param EdgeTol:  (Default value = DBL_TOL)
            EdgeTol:  (Default value = DBL_TOL)

        Returns:

        """

        def _isInsideAndNotOnPosEdge(P, EdgeTol=EdgeTol):
            return (
                self.TileShape.isInShape(P)
                and self.TileShape.satisfiesFacet(P, 2, -EdgeTol)
                and self.TileShape.satisfiesFacet(P, 3, -EdgeTol)
            )

        if _isInsideAndNotOnPosEdge(P, EdgeTol):
            return True  # Do not transform P, return True
        NegCorner = self.TileShape.V[0]
        PosCorner = self.TileShape.V[7]
        Nx = self.TileShape.FacetNorms[3]
        Ny = self.TileShape.FacetNorms[2]
        while np.inner(P - NegCorner, -Nx) > PlanarTiling.DBL_TOL:
            P += self.TileShape.Vx
        while np.inner(P - PosCorner, Nx) > -PlanarTiling.DBL_TOL:
            P -= self.TileShape.Vx
        while np.inner(P - NegCorner, -Ny) > PlanarTiling.DBL_TOL:
            P += self.TileShape.Vy
        while np.inner(P - PosCorner, Ny) > -PlanarTiling.DBL_TOL:
            P -= self.TileShape.Vy
        assert _isInsideAndNotOnPosEdge(P)
        return True

    def getFractionalCoords(
        self, P, blnRelativeToCenter=False, blnRoundInside=True, blnPreferZero=True
    ):
        """

        Args:
            P: param blnRelativeToCenter:  (Default value = False)
            blnRoundInside: Default value = True)
            blnPreferZero: Default value = True)
            blnRelativeToCenter:  (Default value = False)

        Returns:

        """
        # NOTE: The option blnRoundInside is just a way to round numbers close
        #       to the bounds (0,1) so that they never result in values outside
        #       of that range. Some programs (like AtomEye?) are not robust
        #       enough to handle small negative numbers or values greater than 1
        #       in these cases.
        Pfrac = self.TileShape.getFractionalCoords(
            P, blnRelativeToCenter=blnRelativeToCenter
        )
        Pfrac -= Pfrac.astype(int)
        if blnRoundInside:
            if blnPreferZero:
                Pfrac[np.isclose(Pfrac, 0.0, rtol=0.0, atol=PlanarTiling.DBL_TOL)] = 0.0
                Pfrac[np.isclose(Pfrac, 1.0, rtol=0.0, atol=PlanarTiling.DBL_TOL)] = 0.0
            else:
                Pfrac[np.isclose(Pfrac, 0.0, rtol=0.0, atol=PlanarTiling.DBL_TOL)] = 0.0
                Pfrac[np.isclose(Pfrac, 1.0, rtol=0.0, atol=PlanarTiling.DBL_TOL)] = 1.0
        return Pfrac

    def getDistance(self, P0, P1):
        """

        Args:
            P0: param P1:
            P1:

        Returns:

        """
        result = np.linalg.norm(P1 - P0)
        for TilingDirection in self.TilingDirections:
            P1Tiled = P1 + TilingDirection
            TiledDistance = np.linalg.norm(P1Tiled - P0)
            if TiledDistance < result:
                result = TiledDistance
        return result

    def replicateDesign(self, D, nTiles, OldToNewIndices=None, AuxPropMap=None):
        """Create a larger Design by tiling a smaller one.

        Args:
            D (Design): A smaller Design to tile.
            nTiles (int/numpy.ndarray): A specifier for number of
                tiles to replicate in periodic directions.
            OldToNewIndices (dict<int,list<int>>): Optional, a
                dictionary to store the mapping of tiled points
                indices from the smaller Design to the larger.
                (Default value = None)
            AuxPropMap (dict<tuple<string,string>,dict<int,float>>):
                Optional, a mapping of locations to properties.
                NOTE: It is modified in place, so make a copy before.
                (Default value = None)

        Returns:
            (Design): A larger, periodic Design

        """
        if isinstance(nTiles, int):
            nTiles = np.array([nTiles] * 2, dtype=int)
        if OldToNewIndices is None and AuxPropMap is not None:
            OldToNewIndices = {}
        if OldToNewIndices is not None:
            OldToNewIndices.clear()
            for i in range(len(D)):
                OldToNewIndices[i] = []
        result = Design()
        for nx in range(nTiles[0]):
            for ny in range(nTiles[1]):
                Offset = nx * self.Vx + ny * self.Vy
                for i, P in enumerate(D.Canvas.Points):
                    j = len(result)
                    result.add(P + Offset, D.Contents[i])
                    if OldToNewIndices is not None:
                        OldToNewIndices[i].append(j)
        if AuxPropMap is not None:
            for AuxProp in AuxPropMap:
                for i in OldToNewIndices:
                    for j in OldToNewIndices[i]:
                        AuxPropMap[AuxProp][j] = AuxPropMap[AuxProp][i]
        return result

    # === BASIC QUERY METHODS
    @property
    def TileShape(self):
        """ """
        return self._TileShape

    @property
    def TilingDirections(self):
        """ """
        return self._TilingDirections

    @property
    def Vx(self):
        """ """
        return self._TilingDirections[2]

    @property
    def Vy(self):
        """ """
        return self._TilingDirections[3]


class CubicTiling(Tiling):
    """ """

    DBL_TOL = 1e-5

    # === STANDARD CONSTRUCTOR
    def __init__(self, Parallelepiped_):
        super().__init__()
        self._TileShape = Parallelepiped_
        self._TilingDirections = []
        self._TilingDirections.append(-Parallelepiped_.Vz)
        self._TilingDirections.append(-Parallelepiped_.Vx - Parallelepiped_.Vz)
        self._TilingDirections.append(-Parallelepiped_.Vy - Parallelepiped_.Vz)
        self._TilingDirections.append(Parallelepiped_.Vx - Parallelepiped_.Vz)
        self._TilingDirections.append(Parallelepiped_.Vy - Parallelepiped_.Vz)
        self._TilingDirections.append(
            -Parallelepiped_.Vx - Parallelepiped_.Vy - Parallelepiped_.Vz
        )
        self._TilingDirections.append(
            -Parallelepiped_.Vx + Parallelepiped_.Vy - Parallelepiped_.Vz
        )
        self._TilingDirections.append(
            Parallelepiped_.Vx - Parallelepiped_.Vy - Parallelepiped_.Vz
        )
        self._TilingDirections.append(
            Parallelepiped_.Vx + Parallelepiped_.Vy - Parallelepiped_.Vz
        )
        self._TilingDirections.append(-Parallelepiped_.Vx)
        self._TilingDirections.append(-Parallelepiped_.Vy)
        self._TilingDirections.append(Parallelepiped_.Vx)
        self._TilingDirections.append(Parallelepiped_.Vy)
        self._TilingDirections.append(-Parallelepiped_.Vx - Parallelepiped_.Vy)
        self._TilingDirections.append(-Parallelepiped_.Vx + Parallelepiped_.Vy)
        self._TilingDirections.append(Parallelepiped_.Vx - Parallelepiped_.Vy)
        self._TilingDirections.append(Parallelepiped_.Vx + Parallelepiped_.Vy)
        self._TilingDirections.append(Parallelepiped_.Vz)
        self._TilingDirections.append(-Parallelepiped_.Vx + Parallelepiped_.Vz)
        self._TilingDirections.append(-Parallelepiped_.Vy + Parallelepiped_.Vz)
        self._TilingDirections.append(Parallelepiped_.Vx + Parallelepiped_.Vz)
        self._TilingDirections.append(Parallelepiped_.Vy + Parallelepiped_.Vz)
        self._TilingDirections.append(
            -Parallelepiped_.Vx - Parallelepiped_.Vy + Parallelepiped_.Vz
        )
        self._TilingDirections.append(
            -Parallelepiped_.Vx + Parallelepiped_.Vy + Parallelepiped_.Vz
        )
        self._TilingDirections.append(
            Parallelepiped_.Vx - Parallelepiped_.Vy + Parallelepiped_.Vz
        )
        self._TilingDirections.append(
            Parallelepiped_.Vx + Parallelepiped_.Vy + Parallelepiped_.Vz
        )

    # === CONSTRUCTOR - From POSCAR file
    @classmethod
    def fromPOSCAR(cls, filename):
        """

        Args:
            filename:

        Returns:

        """
        return cls(Parallelepiped.fromPOSCAR(filename))

    # === PROPERTY EVALUATION METHODS
    def transformInsideTile(self, P, EdgeTol=DBL_TOL):
        """

        Args:
            P: param EdgeTol:  (Default value = DBL_TOL)
            EdgeTol:  (Default value = DBL_TOL)

        Returns:

        """

        def _isInsideAndNotOnPosEdge(P, EdgeTol=EdgeTol):
            return (
                self.TileShape.isInShape(P)
                and self.TileShape.satisfiesFacet(P, 2, -EdgeTol)
                and self.TileShape.satisfiesFacet(P, 3, -EdgeTol)
                and self.TileShape.satisfiesFacet(P, 4, -EdgeTol)
            )

        if _isInsideAndNotOnPosEdge(P):
            return True  # Do not transform P, return True
        NegCorner = self.TileShape.V[0]
        PosCorner = self.TileShape.V[7]
        Nx = self.TileShape.FacetNorms[3]
        Ny = self.TileShape.FacetNorms[2]
        Nz = self.TileShape.FacetNorms[4]
        while np.inner(P - NegCorner, -Nx) > CubicTiling.DBL_TOL:
            P += self.TileShape.Vx
        while np.inner(P - PosCorner, Nx) > -CubicTiling.DBL_TOL:
            P -= self.TileShape.Vx
        while np.inner(P - NegCorner, -Ny) > CubicTiling.DBL_TOL:
            P += self.TileShape.Vy
        while np.inner(P - PosCorner, Ny) > -CubicTiling.DBL_TOL:
            P -= self.TileShape.Vy
        while np.inner(P - NegCorner, -Nz) > CubicTiling.DBL_TOL:
            P += self.TileShape.Vz
        while np.inner(P - PosCorner, Nz) > -CubicTiling.DBL_TOL:
            P -= self.TileShape.Vz
        assert _isInsideAndNotOnPosEdge(P)
        return True

    def getFractionalCoords(
        self, P, blnRelativeToCenter=False, blnRoundInside=True, blnPreferZero=True
    ):
        """

        Args:
            P: param blnRelativeToCenter:  (Default value = False)
            blnRoundInside: Default value = True)
            blnPreferZero: Default value = True)
            blnRelativeToCenter:  (Default value = False)

        Returns:

        """
        # NOTE: The option blnRoundInside is just a way to round numbers close
        #       to the bounds (0,1) so that they never result in values outside
        #       of that range. Some programs (like AtomEye?) are not robust
        #       enough to handle small negative numbers or values greater than 1
        #       in these cases.
        Pfrac = self.TileShape.getFractionalCoords(
            P, blnRelativeToCenter=blnRelativeToCenter
        )
        Pfrac -= Pfrac.astype(int)
        if blnRoundInside:
            if blnPreferZero:
                Pfrac[np.isclose(Pfrac, 0.0, rtol=0.0, atol=CubicTiling.DBL_TOL)] = 0.0
                Pfrac[np.isclose(Pfrac, 1.0, rtol=0.0, atol=CubicTiling.DBL_TOL)] = 0.0
            else:
                Pfrac[np.isclose(Pfrac, 0.0, rtol=0.0, atol=CubicTiling.DBL_TOL)] = 0.0
                Pfrac[np.isclose(Pfrac, 1.0, rtol=0.0, atol=CubicTiling.DBL_TOL)] = 1.0
        return Pfrac

    def getDistance(self, P0, P1):
        """

        Args:
            P0: param P1:
            P1:

        Returns:

        """
        result = np.linalg.norm(P1 - P0)
        for TilingDirection in self.TilingDirections:
            P1Tiled = P1 + TilingDirection
            TiledDistance = np.linalg.norm(P1Tiled - P0)
            if TiledDistance < result:
                result = TiledDistance
        return result

    def replicateDesign(self, D, nTiles, OldToNewIndices=None, AuxPropMap=None):
        """Create a larger Design by tiling a smaller one.

        Args:
            D (Design): A smaller Design to tile.
            nTiles (int/numpy.ndarray): A specifier for number of
                tiles to replicate in periodic directions.
            OldToNewIndices (dict<int,list<int>>): Optional, a
                dictionary to store the mapping of tiled points
                indices from the smaller Design to the larger.
                (Default value = None)
            AuxPropMap (dict<tuple<string,string>,dict<int,float>>):
                Optional, a mapping of locations to properties.
                NOTE: It is modified in place, so make a copy before.
                (Default value = None)

        Returns:
            (Design): A larger, periodic Design

        """
        if isinstance(nTiles, int):
            nTiles = np.array([nTiles] * 3, dtype=int)
        if OldToNewIndices is None and AuxPropMap is not None:
            OldToNewIndices = {}
        if OldToNewIndices is not None:
            OldToNewIndices.clear()
            for i in range(len(D)):
                OldToNewIndices[i] = []
        result = Design()
        for nx in range(nTiles[0]):
            for ny in range(nTiles[1]):
                for nz in range(nTiles[2]):
                    Offset = nx * self.Vx + ny * self.Vy + nz * self.Vz
                    for i, P in enumerate(D.Canvas.Points):
                        j = len(result)
                        result.add(P + Offset, D.Contents[i])
                        if OldToNewIndices is not None:
                            OldToNewIndices[i].append(j)
        if AuxPropMap is not None:
            for AuxProp in AuxPropMap:
                for i in OldToNewIndices:
                    for j in OldToNewIndices[i]:
                        AuxPropMap[AuxProp][j] = AuxPropMap[AuxProp][i]
        return result

    # === BASIC QUERY METHODS
    @property
    def TileShape(self):
        """ """
        return self._TileShape

    @property
    def TilingDirections(self):
        """ """
        return self._TilingDirections

    @property
    def Vx(self):
        """ """
        return self._TilingDirections[11]

    @property
    def Vy(self):
        """ """
        return self._TilingDirections[12]

    @property
    def Vz(self):
        """ """
        return self._TilingDirections[0]
