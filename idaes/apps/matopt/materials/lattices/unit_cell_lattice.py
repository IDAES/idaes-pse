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
from abc import abstractmethod
from copy import deepcopy
import numpy as np

from ...util.util import myArrayEq
from .lattice import Lattice
from ..transform_func import TransformFunc


class UnitCell(object):
    DBL_TOL = 1e-5

    # === STANDARD CONSTRUCTOR
    def __init__(self, argTiling, FracPositions):
        self._Tiling = argTiling
        self._FracPositions = FracPositions

    # === ASSERTION OF CLASS DESIGN
    def isConsistentWithDesign(self):
        for Position in self.FracPositions:
            for coord in Position:
                if coord < 0.0 - UnitCell.DBL_TOL or coord > 1.0 + UnitCell.DBL_TOL:
                    return False
        return True

    # === MANIPULATION METHODS
    def addPosition(self, P):
        self._FracPositions.append(P)

    def applyTransF(self, TransF):
        if isinstance(TransF, TransformFunc):
            self._Tiling.applyTransF(TransF)
        else:
            raise TypeError

    # === PROPERTY EVALUATION METHODS
    @property
    def Tiling(self):
        return self._Tiling

    @property
    def FracPositions(self):
        return self._FracPositions

    def convertToFrac(self, P, blnDiscardIntPart=True):
        P[:] = self.Tiling.getFractionalCoords(
            P, blnRelativeToCenter=False, blnRoundInside=True, blnPreferZero=True
        )
        if blnDiscardIntPart:
            P -= P.astype(int)

    def getConvertToFrac(self, P, blnDiscardIntPart=True):
        result = deepcopy(P)
        self.convertToFrac(result, blnDiscardIntPart)
        return result

    def getPointType(self, P):
        PFrac = self.getConvertToFrac(P, blnDiscardIntPart=True)
        for i, FracPosition in enumerate(self.FracPositions):
            if myArrayEq(FracPosition, PFrac, UnitCell.DBL_TOL):
                return i
        return None


class UnitCellLattice(Lattice):
    # === AUXILIARY METHODS
    def ScanRef(self, RefScanMin, RefScanMax):
        for n1 in range(RefScanMin[0], RefScanMax[0] + 1):
            for n2 in range(RefScanMin[1], RefScanMax[1] + 1):
                for n3 in range(RefScanMin[2], RefScanMax[2] + 1):
                    for FracPart in self.RefUnitCell.FracPositions:
                        yield (
                            n1 * self.RefUnitCell.Tiling.TileShape.Vx
                            + n2 * self.RefUnitCell.Tiling.TileShape.Vy
                            + n3 * self.RefUnitCell.Tiling.TileShape.Vz
                            + FracPart
                        )

    def Scan(self, argPolyhedron):
        RefScanMin = self._getConvertToReference(argPolyhedron.V[0])
        RefScanMax = self._getConvertToReference(argPolyhedron.V[0])
        for v in argPolyhedron.V:
            RefScanMin = np.minimum(RefScanMin, self._getConvertToReference(v))
            RefScanMax = np.maximum(RefScanMax, self._getConvertToReference(v))
        RefScanMin = RefScanMin.astype(int) + np.array([-1, -1, -1], dtype=int)
        RefScanMax = RefScanMax.astype(int) + np.array([1, 1, 1], dtype=int)
        for RefP in self.ScanRef(RefScanMin, RefScanMax):
            yield self._getConvertFromReference(RefP)

    # === STANDARD CONSTRUCTOR
    def __init__(self, RefUnitCell):
        self._RefUnitCell = RefUnitCell
        Lattice.__init__(self)

    # === PROPERTY EVALUATION METHODS
    def isOnLattice(self, P):
        return (
            self.RefUnitCell.getPointType(Lattice._getConvertToReference(self, P))
            is not None
        )

    @abstractmethod
    def areNeighbors(self, P1, P2):
        raise NotImplementedError

    @abstractmethod
    def getNeighbors(self, P, layer):
        raise NotImplementedError

    # === BASIC QUERY METHODS
    @property
    def RefUnitCell(self):
        return self._RefUnitCell

    @property
    def UnitCell(self):
        result = deepcopy(self.RefUnitCell)
        for f in self._TransformFuncs:
            result.applyTransF(f)
        return result
