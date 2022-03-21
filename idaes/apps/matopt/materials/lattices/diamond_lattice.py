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
from ..geometry import Cube
from ..tiling import CubicTiling
from ..transform_func import ScaleFunc, RotateFunc
from ...util.util import ListHasPoint


class DiamondLattice(UnitCellLattice):
    RefIAD = sqrt(3) / 4

    # === STANDARD CONSTRUCTOR
    def __init__(self, IAD):
        RefUnitCellShape = Cube(1, BotBackLeftCorner=np.array([0, 0, 0], dtype=float))
        RefUnitCellTiling = CubicTiling(RefUnitCellShape)
        RefFracPositions = [
            np.array([0.0, 0.0, 0.0]),
            np.array([0.5, 0.5, 0.0]),
            np.array([0.0, 0.5, 0.5]),
            np.array([0.5, 0.0, 0.5]),
            np.array([0.25, 0.25, 0.25]),
            np.array([0.25, 0.75, 0.75]),
            np.array([0.75, 0.25, 0.75]),
            np.array([0.75, 0.75, 0.25]),
        ]
        RefUnitCell = UnitCell(RefUnitCellTiling, RefFracPositions)
        UnitCellLattice.__init__(self, RefUnitCell)
        self._IAD = (
            DiamondLattice.RefIAD
        )  # IAD is set correctly after calling applyTransF
        self.applyTransF(ScaleFunc(IAD / DiamondLattice.RefIAD))
        self._NthNeighbors = [
            [
                [
                    np.array([0.25, 0.25, 0.25]),
                    np.array([-0.25, -0.25, 0.25]),
                    np.array([-0.25, 0.25, -0.25]),
                    np.array([0.25, -0.25, -0.25]),
                ],
                [
                    np.array([-0.25, -0.25, -0.25]),
                    np.array([0.25, 0.25, -0.25]),
                    np.array([0.25, -0.25, 0.25]),
                    np.array([-0.25, 0.25, 0.25]),
                ],
            ],
            [
                [
                    np.array([0.0, 0.5, 0.5]),
                    np.array([0.0, 0.5, -0.5]),
                    np.array([0.0, -0.5, 0.5]),
                    np.array([0.0, -0.5, -0.5]),
                    np.array([0.5, 0.5, 0.0]),
                    np.array([0.5, 0.0, 0.5]),
                    np.array([0.5, -0.5, 0.0]),
                    np.array([0.5, 0.0, -0.5]),
                    np.array([-0.5, 0.5, 0.0]),
                    np.array([-0.5, 0.0, 0.5]),
                    np.array([-0.5, -0.5, 0.0]),
                    np.array([-0.5, 0.0, -0.5]),
                ],
                [
                    np.array([0.0, 0.5, 0.5]),
                    np.array([0.0, 0.5, -0.5]),
                    np.array([0.0, -0.5, 0.5]),
                    np.array([0.0, -0.5, -0.5]),
                    np.array([0.5, 0.5, 0.0]),
                    np.array([0.5, 0.0, 0.5]),
                    np.array([0.5, -0.5, 0.0]),
                    np.array([0.5, 0.0, -0.5]),
                    np.array([-0.5, 0.5, 0.0]),
                    np.array([-0.5, 0.0, 0.5]),
                    np.array([-0.5, -0.5, 0.0]),
                    np.array([-0.5, 0.0, -0.5]),
                ],
            ],
        ]
        self._typeDict = {0: 0, 3: 1}
        self._relativePositions = {
            0: np.array([0.0, 0.0, 0.0]),
            3: np.array([0.25, 0.25, 0.25]),
        }

    # === CONSTRUCTOR - Aligned with {100}
    @classmethod
    def alignedWith100(cls, IAD):
        return cls(IAD)  # Default implementation

    # === CONSTRUCTOR - Aligned with {110}
    @classmethod
    def aligndWith110(cls, IAD):
        result = cls(IAD)
        thetaX = 0
        thetaY = np.pi * 0.25
        thetaZ = 0
        result.applyTransF(RotateFunc.fromXYZAngles(thetaX, thetaY, thetaZ))
        return result

    # === CONSTRUCTOR - Aligned with {111}
    @classmethod
    def alignedWith111(cls, IAD, blnTrianglesAlignedWithX=True):
        result = cls(IAD)
        thetaX = -np.pi * 0.25
        thetaY = -np.arctan2(-sqrt(2), 2)
        thetaZ = np.pi * 0.5 if blnTrianglesAlignedWithX else 0
        result.applyTransF(RotateFunc.fromXYZAngles(thetaX, thetaY, thetaZ))
        return result

    # === CONSTRUCTOR - Aligned with {xyz}
    @classmethod
    def alignedWith(cls, IAD, MI):
        if (type(MI) is str) and (len(MI) == 3) and all(x.isdigit() for x in MI):
            if MI in ["100", "010", "001"]:
                return cls(IAD)
            elif MI in ["110", "101", "011"]:
                return cls.aligndWith110(IAD)
            elif MI == "111":
                return cls.alignedWith111(IAD)
            else:
                result = cls(IAD)
                a = np.array([0.0, 0.0, 1.0])
                b = np.array([float(MI[0]), float(MI[1]), float(MI[2])])
                axis = np.cross(a, b)
                angle = np.arccos(
                    np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))
                )
                result.applyTransF(RotateFunc.fromAxisAngle(axis, angle))
                return result
        return ValueError("DiamondLattice.alignedWith: Input direction is not correct.")

    # === MANIPULATION METHODS
    def applyTransF(self, TransF):
        if isinstance(TransF, ScaleFunc):
            if TransF.isIsometric:
                self._IAD *= TransF.Scale[0]
            else:
                raise ValueError(
                    "DiamondLattice.applyTransF: Can only scale isometrically"
                )
        UnitCellLattice.applyTransF(self, TransF)

    # === AUXILIARY METHODS
    def _getPointType(self, P):
        return (int(round(P[0] * 4)) + int(round(P[1] * 4)) + int(round(P[2] * 4))) % 4

    # === PROPERTY EVALUATION METHODS
    # NOTE: inherited from UnitCellLattice
    # def isOnLattice(self,P):

    def areNeighbors(self, P1, P2):
        return np.linalg.norm(P2 - P1) <= self.IAD

    def getNeighbors(self, P, layer=1):
        RefP = self._getConvertToReference(P)
        PType = self._getPointType(RefP)
        if PType not in self._typeDict.keys():
            raise ValueError("DiamondLattice.getNeighbors Should never reach here!")
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
                        if not ListHasPoint(NList[v], N, 0.001 * DiamondLattice.RefIAD):
                            tmp[v].append(N)
                            NList[v].append(N)
            self._NthNeighbors.append(tmp)

    def isASite(self, P):
        RefP = self._getConvertToReference(P)
        PType = self._getPointType(RefP)
        return PType == 0

    def isBSite(self, P):
        RefP = self._getConvertToReference(P)
        PType = self._getPointType(RefP)
        return PType == 3

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

    @property
    def Diamond100LayerSpacing(self):
        return self.IAD / sqrt(3)

    @property
    def Diamond110LayerSpacing(self):
        return self.IAD * sqrt(2) / sqrt(3)

    @property
    def Diamond111LayerSpacing(self):
        return self.IAD * 4 / 3

    @property
    def Diamond112LayerSpacing(self):
        return self.IAD * sqrt(2) / 3

    def getLayerSpacing(self, MI):
        if (type(MI) is str) and (len(MI) == 3) and all(x.isdigit() for x in MI):
            if MI in ["100", "010", "001"]:
                return self.Diamond100LayerSpacing
            elif MI in ["110", "101", "011"]:
                return self.Diamond110LayerSpacing
            elif MI == "111":
                return self.Diamond111LayerSpacing
            elif MI in ["112", "121", "211"]:
                return self.Diamond112LayerSpacing
            else:
                raise NotImplementedError(
                    "DiamondLattice.getLayerSpacing: Input direction is not supported."
                )
        return ValueError(
            "DiamondLattice.getLayerSpacing: Input direction is not correct."
        )

    def getShellSpacing(self, MI):
        if (type(MI) is str) and (len(MI) == 3) and all(x.isdigit() for x in MI):
            if MI in ["100", "010", "001", "110", "101", "011", "111"]:
                return self.IAD * sqrt(8) / sqrt(3)
            elif MI in ["112", "121", "211"]:
                return self.IAD * sqrt(2) / sqrt(3)
            else:
                raise NotImplementedError(
                    "DiamondLattice.getShellSpacing: Input direction is not supported."
                )
        return ValueError("The input direction is not correct.")

    def getUniqueLayerCount(self, MI):
        if (type(MI) is str) and (len(MI) == 3) and all(x.isdigit() for x in MI):
            if MI in ["100", "010", "001"]:
                return 4
            elif MI in ["110", "101", "011"]:
                return 2
            elif MI == "111":
                return 3
            elif MI in ["112", "121", "211"]:
                return 6
            else:
                raise NotImplementedError(
                    "DiamondLattice.getUniqueLayerCount: Input direction is not supported."
                )
        return ValueError("The input direction is not correct.")
