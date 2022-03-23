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

import numpy as np

from .unit_cell_lattice import UnitCell, UnitCellLattice
from ..geometry import RectPrism
from ..tiling import CubicTiling
from ..transform_func import ScaleFunc, RotateFunc, ReflectFunc


class PerovskiteLattice(UnitCellLattice):
    RefA = 1
    RefB = 1
    RefC = 1

    # === STANDARD CONSTRUCTOR
    def __init__(self, A, B, C):
        RefUnitCellShape = RectPrism(
            PerovskiteLattice.RefA,
            PerovskiteLattice.RefB,
            PerovskiteLattice.RefC,
            np.array([0, 0, 0], dtype=float),
        )
        RefUnitCellTiling = CubicTiling(RefUnitCellShape)
        RefFracPositions = [
            np.array([0.0, 0.0, 0.0]),
            np.array([0.5, 0.5, 0.5]),
            np.array([0.5, 0.5, 0.0]),
            np.array([0.5, 0.0, 0.5]),
            np.array([0.0, 0.5, 0.5]),
        ]
        RefUnitCell = UnitCell(RefUnitCellTiling, RefFracPositions)
        UnitCellLattice.__init__(self, RefUnitCell)
        self._A = PerovskiteLattice.RefA
        self._B = PerovskiteLattice.RefB
        self._C = PerovskiteLattice.RefC
        self.applyTransF(
            ScaleFunc(
                np.array(
                    [
                        A / PerovskiteLattice.RefA,
                        B / PerovskiteLattice.RefB,
                        C / PerovskiteLattice.RefC,
                    ]
                )
            )
        )

    # === MANIPULATION METHODS
    def applyTransF(self, TransF):
        if isinstance(TransF, ScaleFunc):
            self._A *= TransF.Scale[0]
            self._B *= TransF.Scale[1]
            self._C *= TransF.Scale[2]
        UnitCellLattice.applyTransF(self, TransF)

    # === PROPERTY EVALUATION METHODS
    # def isOnLattice(self,P):

    def areNeighbors(self, P1, P2):
        raise NotImplementedError(
            "PerovskiteLattice: there is no universal definition of nearest neighbors."
        )

    def getNeighbors(self, P, layer=1):
        if layer > 2:
            raise ValueError(
                "PerovskiteLattice: there is no universal definition of N-th nearest neighbors."
            )
        RefP = self._getConvertToReference(P)
        if self.isASite(P):
            return []
        elif self.isBSite(P):
            RefNeighs = [
                np.array([0.5, 0.0, 0.0]),
                np.array([0.0, 0.5, 0.0]),
                np.array([0.0, 0.0, 0.5]),
                np.array([-0.5, 0.0, 0.0]),
                np.array([0.0, -0.5, 0.0]),
                np.array([0.0, 0.0, -0.5]),
            ]
        elif self.isOSite(P):
            PointType = self.RefUnitCell.getPointType(self._getConvertToReference(P))
            if PointType == 4:  # i.e., motif aligned with x-axis
                RefNeighs = [
                    np.array([-0.5, 0.0, 0.0]),
                    np.array([0.5, 0.0, 0.0]),
                    np.array([-0.5, 1.0, 0.0]),
                    np.array([-0.5, 0.0, 1.0]),
                    np.array([-0.5, -1.0, 0.0]),
                    np.array([-0.5, 0.0, -1.0]),
                    np.array([0.5, 1.0, 0.0]),
                    np.array([0.5, 0.0, 1.0]),
                    np.array([0.5, -1.0, 0.0]),
                    np.array([0.5, 0.0, -1.0]),
                ]
            elif PointType == 3:  # i.e., motif aligned with y-axis
                RefNeighs = [
                    np.array([0.0, -0.5, 0.0]),
                    np.array([0.0, 0.5, 0.0]),
                    np.array([1.0, -0.5, 0.0]),
                    np.array([0.0, -0.5, 1.0]),
                    np.array([-1.0, -0.5, 0.0]),
                    np.array([0.0, -0.5, -1.0]),
                    np.array([1.0, 0.5, 0.0]),
                    np.array([0.0, 0.5, 1.0]),
                    np.array([-1.0, 0.5, 0.0]),
                    np.array([0.0, 0.5, -1.0]),
                ]
            elif PointType == 2:  # i.e., motif aligned with z-axis
                RefNeighs = [
                    np.array([0.0, 0.0, -0.5]),
                    np.array([0.0, 0.0, 0.5]),
                    np.array([1.0, 0.0, -0.5]),
                    np.array([0.0, 1.0, -0.5]),
                    np.array([-1.0, 0.0, -0.5]),
                    np.array([0.0, -1.0, -0.5]),
                    np.array([1.0, 0.0, 0.5]),
                    np.array([0.0, 1.0, 0.5]),
                    np.array([-1.0, 0.0, 0.5]),
                    np.array([0.0, -1.0, 0.5]),
                ]
            else:
                raise ValueError(
                    "PerovskiteLattice.getNeighbors Should never reach here!"
                )
        else:
            raise ValueError(
                "PerovskiteLattice.getNeighbors given point apparently not on lattice"
            )
        result = deepcopy(RefNeighs)
        for NeighP in result:
            NeighP += RefP
            self._convertFromReference(NeighP)
        return result

    def isASite(self, P):
        return self.RefUnitCell.getPointType(self._getConvertToReference(P)) == 0

    def isBSite(self, P):
        return self.RefUnitCell.getPointType(self._getConvertToReference(P)) == 1

    def isOSite(self, P):
        PointType = self.RefUnitCell.getPointType(self._getConvertToReference(P))
        return PointType == 2 or PointType == 3 or PointType == 4

    def setDesign(self, D, AType, BType, OType):
        for i, P in enumerate(D.Canvas.Points):
            if self.isASite(P):
                D.setContent(i, AType)
            elif self.isBSite(P):
                D.setContent(i, BType)
            elif self.isOSite(P):
                D.setContent(i, OType)
            else:
                raise ValueError("setDesign can not set site not on lattice")

    # === BASIC QUERY METHODS
    @property
    def A(self):
        return self._A

    @property
    def B(self):
        return self._B

    @property
    def C(self):
        return self._C


def getOxygenSymTransFs():
    result = []
    ReflX = ReflectFunc.acrossX()
    ReflY = ReflectFunc.acrossY()
    ReflZ = ReflectFunc.acrossZ()
    RotYY = RotateFunc.fromXYZAngles(0, np.pi, 0)
    RotZZ = RotateFunc.fromXYZAngles(0, 0, np.pi)
    RotX = RotateFunc.fromXYZAngles(np.pi * 0.5, 0, 0)
    RotXX = RotateFunc.fromXYZAngles(np.pi, 0, 0)
    RotXXX = RotateFunc.fromXYZAngles(np.pi * 1.5, 0, 0)

    result.append(RotX)
    result.append(RotXX)
    result.append(RotXXX)

    result.append(ReflX)
    result.append(ReflX + RotX)
    result.append(ReflX + RotXX)
    result.append(ReflX + RotXXX)

    result.append(ReflY)
    result.append(ReflY + RotX)
    result.append(ReflY + RotXX)
    result.append(ReflY + RotXXX)

    result.append(ReflZ)
    result.append(ReflZ + RotX)
    result.append(ReflZ + RotXX)
    result.append(ReflZ + RotXXX)

    result.append(ReflX + ReflY)
    result.append(ReflX + ReflY + RotX)
    result.append(ReflX + ReflY + RotXX)
    result.append(ReflX + ReflY + RotXXX)

    result.append(ReflX + ReflZ)
    result.append(ReflX + ReflZ + RotX)
    result.append(ReflX + ReflZ + RotXX)
    result.append(ReflX + ReflZ + RotXXX)

    result.append(ReflY + ReflZ)
    result.append(ReflY + ReflZ + RotX)
    result.append(ReflY + ReflZ + RotXX)
    result.append(ReflY + ReflZ + RotXXX)

    result.append(ReflX + ReflY + ReflZ)
    result.append(ReflX + ReflY + ReflZ + RotX)
    result.append(ReflX + ReflY + ReflZ + RotXX)
    result.append(ReflX + ReflY + ReflZ + RotXXX)

    result.append(RotYY + RotX)
    result.append(RotYY + RotXX)
    result.append(RotYY + RotXXX)

    result.append(RotYY + ReflX)
    result.append(RotYY + ReflX + RotX)
    result.append(RotYY + ReflX + RotXX)
    result.append(RotYY + ReflX + RotXXX)

    result.append(RotYY + ReflY)
    result.append(RotYY + ReflY + RotX)
    result.append(RotYY + ReflY + RotXX)
    result.append(RotYY + ReflY + RotXXX)

    result.append(RotYY + ReflZ)
    result.append(RotYY + ReflZ + RotX)
    result.append(RotYY + ReflZ + RotXX)
    result.append(RotYY + ReflZ + RotXXX)

    result.append(RotYY + ReflX + ReflY)
    result.append(RotYY + ReflX + ReflY + RotX)
    result.append(RotYY + ReflX + ReflY + RotXX)
    result.append(RotYY + ReflX + ReflY + RotXXX)

    result.append(RotYY + ReflX + ReflZ)
    result.append(RotYY + ReflX + ReflZ + RotX)
    result.append(RotYY + ReflX + ReflZ + RotXX)
    result.append(RotYY + ReflX + ReflZ + RotXXX)

    result.append(RotYY + ReflY + ReflZ)
    result.append(RotYY + ReflY + ReflZ + RotX)
    result.append(RotYY + ReflY + ReflZ + RotXX)
    result.append(RotYY + ReflY + ReflZ + RotXXX)

    result.append(RotYY + ReflX + ReflY + ReflZ)
    result.append(RotYY + ReflX + ReflY + ReflZ + RotX)
    result.append(RotYY + ReflX + ReflY + ReflZ + RotXX)
    result.append(RotYY + ReflX + ReflY + ReflZ + RotXXX)

    result.append(RotZZ + RotX)
    result.append(RotZZ + RotXX)
    result.append(RotZZ + RotXXX)

    result.append(RotZZ + ReflX)
    result.append(RotZZ + ReflX + RotX)
    result.append(RotZZ + ReflX + RotXX)
    result.append(RotZZ + ReflX + RotXXX)

    result.append(RotZZ + ReflY)
    result.append(RotZZ + ReflY + RotX)
    result.append(RotZZ + ReflY + RotXX)
    result.append(RotZZ + ReflY + RotXXX)

    result.append(RotZZ + ReflZ)
    result.append(RotZZ + ReflZ + RotX)
    result.append(RotZZ + ReflZ + RotXX)
    result.append(RotZZ + ReflZ + RotXXX)

    result.append(RotZZ + ReflX + ReflY)
    result.append(RotZZ + ReflX + ReflY + RotX)
    result.append(RotZZ + ReflX + ReflY + RotXX)
    result.append(RotZZ + ReflX + ReflY + RotXXX)

    result.append(RotZZ + ReflX + ReflZ)
    result.append(RotZZ + ReflX + ReflZ + RotX)
    result.append(RotZZ + ReflX + ReflZ + RotXX)
    result.append(RotZZ + ReflX + ReflZ + RotXXX)

    result.append(RotZZ + ReflY + ReflZ)
    result.append(RotZZ + ReflY + ReflZ + RotX)
    result.append(RotZZ + ReflY + ReflZ + RotXX)
    result.append(RotZZ + ReflY + ReflZ + RotXXX)

    result.append(RotZZ + ReflX + ReflY + ReflZ)
    result.append(RotZZ + ReflX + ReflY + ReflZ + RotX)
    result.append(RotZZ + ReflX + ReflY + ReflZ + RotXX)
    result.append(RotZZ + ReflX + ReflY + ReflZ + RotXXX)

    return result
