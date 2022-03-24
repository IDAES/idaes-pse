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
from copy import deepcopy
from math import sqrt

import numpy as np

from .unit_cell_lattice import UnitCell, UnitCellLattice
from ..geometry import Parallelepiped
from ..tiling import CubicTiling
from ..transform_func import ScaleFunc, RotateFunc
from ...opt import DBL_TOL
from ...util.util import ListHasPoint


class WurtziteLattice(UnitCellLattice):
    RefIAD = sqrt(3 / 8)

    # === STANDARD CONSTRUCTOR
    def __init__(self, IAD):
        RefUnitCellShape = Parallelepiped(
            np.array([1, 0, 0], dtype=float),
            np.array([0.5, sqrt(3) / 2, 0], dtype=float),
            np.array([0, 0, sqrt(8) / sqrt(3)], dtype=float),
            BotBackLeftCorner=np.array([0, 0, 0], dtype=float),
        )
        RefUnitCellTiling = CubicTiling(RefUnitCellShape)
        RefFracPositions = [
            np.array([0.0, 0.0, 0.0]),
            np.array([0.5, sqrt(3) / 6.0, 1.0 / sqrt(24)]),
            np.array([0.5, sqrt(3) / 6.0, 4.0 / sqrt(24)]),
            np.array([0.0, 0.0, 5.0 / sqrt(24)]),
        ]
        RefUnitCell = UnitCell(RefUnitCellTiling, RefFracPositions)
        UnitCellLattice.__init__(self, RefUnitCell)
        self._IAD = (
            WurtziteLattice.RefIAD
        )  # IAD is set correctly after calling applyTransF
        self.applyTransF(ScaleFunc(IAD / WurtziteLattice.RefIAD))
        self._NthNeighbors = [
            [
                [
                    np.array([0.0, 0.0, -3 / sqrt(24)]),
                    np.array([0.0, -sqrt(3) / 3, 1 / sqrt(24)]),
                    np.array([0.5, sqrt(3) / 6, 1 / sqrt(24)]),
                    np.array([-0.5, sqrt(3) / 6, 1 / sqrt(24)]),
                ],
                [
                    np.array([0.0, 0.0, 3 / sqrt(24)]),
                    np.array([0.0, sqrt(3) / 3, -1 / sqrt(24)]),
                    np.array([0.5, -sqrt(3) / 6, -1 / sqrt(24)]),
                    np.array([-0.5, -sqrt(3) / 6, -1 / sqrt(24)]),
                ],
                [
                    np.array([0.0, 0.0, -3 / sqrt(24)]),
                    np.array([0.0, sqrt(3) / 3, 1 / sqrt(24)]),
                    np.array([0.5, -sqrt(3) / 6, 1 / sqrt(24)]),
                    np.array([-0.5, -sqrt(3) / 6, 1 / sqrt(24)]),
                ],
                [
                    np.array([0.0, 0.0, 3 / sqrt(24)]),
                    np.array([0.0, -sqrt(3) / 3, -1 / sqrt(24)]),
                    np.array([0.5, sqrt(3) / 6, -1 / sqrt(24)]),
                    np.array([-0.5, sqrt(3) / 6, -1 / sqrt(24)]),
                ],
            ]
        ]
        self._typeDict = {0: 0, 1: 1, 4: 2, 5: 3}
        self._relativePositions = {
            0: np.array([0.0, 0.0, 0.0]),
            1: np.array([0.5, sqrt(3) / 6.0, 1.0 / sqrt(24)]),
            4: np.array([0.5, sqrt(3) / 6.0, 4.0 / sqrt(24)]),
            5: np.array([0.0, 0.0, 5.0 / sqrt(24)]),
        }

    @classmethod
    def alignedWith(cls, IAD, MI):
        if MI == "0001":
            return cls(IAD)
        elif MI == "1100":
            result = cls(IAD)
            axis = np.array([1, 0, 0], dtype=float)
            angle = np.pi / 2
            result.applyTransF(RotateFunc.fromAxisAngle(axis, angle))
            return result
        elif MI == "1120":
            result = cls(IAD)
            axis = np.array([0, 1, 0], dtype=float)
            angle = np.pi / 2
            result.applyTransF(RotateFunc.fromAxisAngle(axis, angle))
            return result
        else:
            raise ValueError(
                "WurtziteLattice.alignedWith: Input direction is not supported."
            )

    # === MANIPULATION METHODS
    def applyTransF(self, TransF):
        if isinstance(TransF, ScaleFunc):
            if TransF.isIsometric:
                self._IAD *= TransF.Scale[0]
            else:
                raise ValueError(
                    "WurtziteLattice.applyTransF: Can only scale isometrically"
                )
        UnitCellLattice.applyTransF(self, TransF)

    # === AUXILIARY METHODS
    def _getPointType(self, P):
        return int(round(P[2] * sqrt(24))) % 8

    # === PROPERTY EVALUATION METHODS
    # NOTE: inherited from UnitCellLattice
    # def isOnLattice(self,P):

    def areNeighbors(self, P1, P2):
        return np.linalg.norm(P2 - P1) <= self.IAD + DBL_TOL

    def getNeighbors(self, P, layer=1):
        RefP = self._getConvertToReference(P)
        PType = self._getPointType(RefP)
        if PType not in self._typeDict.keys():
            raise ValueError("WurtziteLattice.getNeighbors Should never reach here!")
        if layer > len(self._NthNeighbors):
            self._calculateNeighbors(layer)
        NBs = deepcopy(self._NthNeighbors[layer - 1][self._typeDict[PType]])
        for NeighP in NBs:
            NeighP += RefP
            self._convertFromReference(NeighP)
        return NBs

    def _calculateNeighbors(self, layer):
        NList = []
        for k, v in self._typeDict.items():
            tmp = [np.array([0, 0, 0], dtype=float)]
            for nb in self._NthNeighbors:
                tmp.extend(nb[v])
            NList.append(tmp)
        for _ in range(layer - len(self._NthNeighbors)):
            tmp = [[] for _ in self._typeDict.keys()]
            for k, v in self._typeDict.items():
                for P in self._NthNeighbors[len(self._NthNeighbors) - 1][v]:
                    PType = self._getPointType(P + self._relativePositions[k])
                    for Q in self._NthNeighbors[0][self._typeDict[PType]]:
                        N = P + Q
                        if not ListHasPoint(
                            NList[v], N, 0.001 * WurtziteLattice.RefIAD
                        ):
                            tmp[v].append(N)
                            NList[v].append(N)
            self._NthNeighbors.append(tmp)

    def isASite(self, P):
        RefP = self._getConvertToReference(P)
        PType = self._getPointType(RefP)
        return PType == 0 or PType == 4

    def isBSite(self, P):
        RefP = self._getConvertToReference(P)
        PType = self._getPointType(RefP)
        return PType == 1 or PType == 5

    def setDesign(self, D, AType, BType):
        for i, P in enumerate(D.Canvas.Points):
            if self.isASite(P):
                D.setContent(i, AType)
            elif self.isBSite(P):
                D.setContent(i, BType)
            else:
                raise ValueError("setDesign can not set site not on lattice")

    # === BASIC QUERY METHODS
    @property
    def IAD(self):
        return self._IAD

    def getLayerSpacing(self, MI):
        if MI == "0001":
            return self._IAD * 4 / 3
        elif MI == "1100":
            return self._IAD * sqrt(2)
        elif MI == "1120":
            return self._IAD * sqrt(8 / 3)
        else:
            raise NotImplementedError(
                "WurtziteLattice.getLayerSpacing: Input direction is not supported."
            )

    def getShellSpacing(self, MI):
        if MI == "0001":
            return self._IAD * sqrt(2)
        elif MI == "1100":
            return self._IAD * sqrt(8 / 3)
        elif MI == "1120":
            return self._IAD * sqrt(8 / 3)
        else:
            raise NotImplementedError(
                "WurtziteLattice.getShellSpacing: Input direction is not supported."
            )

    def getUniqueLayerCount(self, MI):
        if MI == "0001":
            return 2
        elif MI == "1100":
            return 2
        elif MI == "1120":
            return 1
        else:
            raise NotImplementedError(
                "WurtziteLattice.getUniqueLayerCount: Input direction is not supported."
            )
