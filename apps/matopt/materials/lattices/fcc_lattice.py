import numpy as np
from math import sqrt
from copy import deepcopy

from ..geometry import Cube, Cuboctahedron
from ..transform_func import ScaleFunc, RotateFunc, ReflectFunc
from .unit_cell_lattice import UnitCell, UnitCellLattice
from ..tiling import CubicTiling


class FCCLattice(UnitCellLattice):
    RefIAD = sqrt(2) / 2

    # === STANDARD CONSTRUCTOR
    def __init__(self, IAD):
        RefUnitCellShape = Cube(1, BotBackLeftCorner=np.array([0, 0, 0], dtype=float))
        RefUnitCellTiling = CubicTiling(RefUnitCellShape)
        RefFracPositions = [np.array([0.0, 0.0, 0.0]),
                            np.array([0.5, 0.5, 0.0]),
                            np.array([0.0, 0.5, 0.5]),
                            np.array([0.5, 0.0, 0.5])]
        RefUnitCell = UnitCell(RefUnitCellTiling, RefFracPositions)
        UnitCellLattice.__init__(self, RefUnitCell)
        self._IAD = FCCLattice.RefIAD  # IAD is set correctly after calling applyTransF
        self._RefNeighborsPattern = [np.array([0.0, -0.5, 0.5]),
                                     np.array([-0.5, -0.5, 0.0]),
                                     np.array([-0.5, 0.0, 0.5]),
                                     np.array([0.5, -0.5, 0.0]),
                                     np.array([0.0, -0.5, -0.5]),
                                     np.array([-0.5, 0.0, -0.5]),
                                     np.array([-0.5, 0.5, 0.0]),
                                     np.array([0.0, 0.5, 0.5]),
                                     np.array([0.5, 0.0, 0.5]),
                                     np.array([0.5, 0.0, -0.5]),
                                     np.array([0.0, 0.5, -0.5]),
                                     np.array([0.5, 0.5, 0.0])]
        self.applyTransF(ScaleFunc(IAD / FCCLattice.RefIAD))
        # import code; code.interact(local=dict(locals(),**globals()));
        assert (self.isConsistentWithDesign())

    # === CONSTRUCTOR - Aligned with FCC {100}
    @classmethod
    def alignedWith100(cls, IAD):
        return cls(IAD)  # Default implementation

    # === CONSTRUCTOR - Aligned with FCC {110}
    @classmethod
    def alignedWith110(cls, IAD):
        result = cls(IAD)
        raise NotImplementedError('TODO')
        return result

    # === CONSTRUCTOR - Aligned with FCC {111}
    @classmethod
    def alignedWith111(cls, IAD, blnTrianglesAlignedWithX=True):
        result = cls(IAD)
        thetaX = -np.pi * 0.25
        thetaY = -np.arctan2(-sqrt(2), 2)
        thetaZ = (np.pi * 0.5 if blnTrianglesAlignedWithX else 0)
        result.applyTransF(RotateFunc.fromXYZAngles(thetaX, thetaY, thetaZ))
        return result

    # === ASSERTION OF CLASS DESIGN
    def isConsistentWithDesign(self):
        # TODO: Put some real checks here
        return UnitCellLattice.isConsistentWithDesign(self)

    # === MANIPULATION METHODS
    def applyTransF(self, TransF):
        if (isinstance(TransF, ScaleFunc)):
            if (TransF.isIsometric):
                self._IAD *= TransF.Scale[0]
            else:
                raise ValueError('FCCLattice applyTransF: Can only scale isometrically')
        UnitCellLattice.applyTransF(self, TransF)

    # === PROPERTY EVALUATION METHODS
    # def isOnLattice(self,P):

    def areNeighbors(self, P1, P2):
        # return distEucPoint2Point(P1,P2) <= self.IAD
        return np.linalg.norm(P2 - P1) <= self.IAD

    def getNeighbors(self, P):
        RefP = self._getConvertToReference(P)
        result = deepcopy(self._RefNeighborsPattern)
        # print(result)
        for NeighP in result:
            NeighP += RefP
            self._convertFromReference(NeighP)
        # print
        # print(result)
        # import code; code.interact(local=dict(locals(),**globals()));
        return result

    # === BASIC QUERY METHODS
    @property
    def IAD(self):
        return self._IAD

    @property
    def FCC111LayerSpacing(self):
        return self.IAD * sqrt(2) / sqrt(3)

    @property
    def FCC100LayerSpacing(self):
        return self.IAD * 0.5

    @property
    def FCC110LayerSpacing(self):
        return self.IAD * sqrt(2) / 2


def getPossibleRotAndRefls111():
    result = []
    thetaX = -np.pi * 0.25
    thetaY = -np.arctan2(-sqrt(2), 2)
    thetaZ = np.pi * 0.5
    S = Cuboctahedron(1)
    ThreeFoldAxis1 = S.FacetDirections[6]
    ThreeFoldAxis2 = S.FacetDirections[7]
    ThreeFoldAxis3 = S.FacetDirections[8]
    ThreeFoldAxis4 = S.FacetDirections[9]
    FourFoldAxis1 = S.FacetDirections[0]
    FourFoldAxis2 = S.FacetDirections[1]
    FourFoldAxis3 = S.FacetDirections[2]

    Rot1ThreeFold1 = RotateFunc.fromAxisAngle(ThreeFoldAxis1, np.pi * 2.0 / 3.0)
    Rot2ThreeFold1 = RotateFunc.fromAxisAngle(ThreeFoldAxis1, np.pi * 4.0 / 3.0)
    Rot1ThreeFold2 = RotateFunc.fromAxisAngle(ThreeFoldAxis2, np.pi * 2.0 / 3.0)
    Rot2ThreeFold2 = RotateFunc.fromAxisAngle(ThreeFoldAxis2, np.pi * 4.0 / 3.0)
    Rot1ThreeFold3 = RotateFunc.fromAxisAngle(ThreeFoldAxis3, np.pi * 2.0 / 3.0)
    Rot2ThreeFold3 = RotateFunc.fromAxisAngle(ThreeFoldAxis3, np.pi * 4.0 / 3.0)
    Rot1ThreeFold4 = RotateFunc.fromAxisAngle(ThreeFoldAxis4, np.pi * 2.0 / 3.0)
    Rot2ThreeFold4 = RotateFunc.fromAxisAngle(ThreeFoldAxis4, np.pi * 4.0 / 3.0)

    Rot1FourFold1 = RotateFunc.fromAxisAngle(FourFoldAxis1, np.pi * 1.0 / 2.0)
    Rot2FourFold1 = RotateFunc.fromAxisAngle(FourFoldAxis1, np.pi * 2.0 / 2.0)
    Rot3FourFold1 = RotateFunc.fromAxisAngle(FourFoldAxis1, np.pi * 3.0 / 2.0)
    Rot1FourFold2 = RotateFunc.fromAxisAngle(FourFoldAxis2, np.pi * 1.0 / 2.0)
    Rot2FourFold2 = RotateFunc.fromAxisAngle(FourFoldAxis2, np.pi * 2.0 / 2.0)
    Rot3FourFold2 = RotateFunc.fromAxisAngle(FourFoldAxis2, np.pi * 3.0 / 2.0)
    Rot1FourFold3 = RotateFunc.fromAxisAngle(FourFoldAxis3, np.pi * 1.0 / 2.0)
    Rot2FourFold3 = RotateFunc.fromAxisAngle(FourFoldAxis3, np.pi * 2.0 / 2.0)
    Rot3FourFold3 = RotateFunc.fromAxisAngle(FourFoldAxis3, np.pi * 3.0 / 2.0)

    ReflX = ReflectFunc.acrossX()
    ReflY = ReflectFunc.acrossY()
    ReflZ = ReflectFunc.acrossZ()
    ReflXY = ReflectFunc(np.array([sqrt(2), sqrt(2), 0]), np.array([0, 0, 0], dtype=float))
    ReflXZ = ReflectFunc(np.array([sqrt(2), 0, sqrt(2)]), np.array([0, 0, 0], dtype=float))
    ReflYZ = ReflectFunc(np.array([0, sqrt(2), sqrt(2)]), np.array([0, 0, 0], dtype=float))
    ReflXYneg = ReflectFunc(np.array([-sqrt(2), sqrt(2), 0]), np.array([0, 0, 0], dtype=float))
    ReflXZneg = ReflectFunc(np.array([-sqrt(2), 0, sqrt(2)]), np.array([0, 0, 0], dtype=float))
    ReflYZneg = ReflectFunc(np.array([0, -sqrt(2), sqrt(2)]), np.array([0, 0, 0], dtype=float))

    raise NotImplementedError('TODO')
