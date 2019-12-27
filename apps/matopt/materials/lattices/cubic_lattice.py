import numpy as np
from math import sqrt
from copy import deepcopy

from ..geometry import Parallelepiped
from ..transform_func import ScaleFunc
from .unit_cell_lattice import UnitCell, UnitCellLattice
from ..tiling import CubicTiling


class CubicLattice(UnitCellLattice):
    RefIAD = 1

    '''
    # === AUXILIARY METHODS
    def __isValidScanPoint(self,P):
        return True # all integer points are on the cubic lattice
    '''

    '''
    def _setupScan(self,ScanMin,ScanMax):
        self.__ScanMin = ScanMin.astype(int)+np.array([-1,-1,-1],dtype=int)
        self.__ScanMax = ScanMax.astype(int)+np.array([ 1, 1, 1],dtype=int)
        self.__ScanP = deepcopy(self.__ScanMin)
        self.__ScanP[0] -= 1
        #print('Setting up ScanMin={}'.format(self.__ScanMin))
        #print('Setting up ScanMax={}'.format(self.__ScanMax))
        #print('Setting up ScanP={}'.format(self.__ScanP))
    def _getNextPoint(self):
        while(self.__ScanP[2]<=self.__ScanMax[2]):
            while(self.__ScanP[1]<=self.__ScanMax[1]):
                while(self.__ScanP[0]<=self.__ScanMax[0]):
                    self.__ScanP[0] += 1
                    if(self.__ScanP[0]<=self.__ScanMax[0] and
                       self.__ScanP[1]<=self.__ScanMax[1] and
                       self.__ScanP[2]<=self.__ScanMax[2]):
                        #print('Returning P={}'.format(self._getConvertFromReference(self.__ScanP.astype(float))))
                        return self._getConvertFromReference(self.__ScanP.astype(float))
                self.__ScanP[0] = self.__ScanMin[0]
                self.__ScanP[1] += 1
            self.__ScanP[1] = self.__ScanMin[1]
            self.__ScanP[2] += 1
        self.__ScanP[2] = self.__ScanMin[2]
        return None # only invoked when ScanP is reset to ScanMin
    '''

    # === STANDARD CONSTRUCTOR
    def __init__(self, IAD):
        RefUnitCellShape = Parallelepiped.fromEdgesAndAngles(CubicLattice.RefIAD,
                                                             CubicLattice.RefIAD,
                                                             CubicLattice.RefIAD,
                                                             np.pi / 2, np.pi / 2, np.pi / 2,
                                                             np.array([0, 0, 0], dtype=float))
        RefUnitCellTiling = CubicTiling(RefUnitCellShape)
        RefFracPositions = [np.array([0.0, 0.0, 0.0])]
        RefUnitCell = UnitCell(RefUnitCellTiling, RefFracPositions)
        UnitCellLattice.__init__(self, RefUnitCell)
        self._IAD = IAD
        self._RefNeighborsPattern = [np.array([1.0, 0.0, 0.0]),
                                     np.array([0.0, 1.0, 0.0]),
                                     np.array([0.0, 0.0, 1.0]),
                                     np.array([-1.0, 0.0, 0.0]),
                                     np.array([0.0, -1.0, 0.0]),
                                     np.array([0.0, 0.0, -1.0])]
        self.applyTransF(ScaleFunc(IAD / CubicLattice.RefIAD))
        assert (self.isConsistentWithDesign())

    # === ASSERTION OF CLASS DESIGN
    def isConsistentWithDesign(self):
        # TODO: Put some real checks here
        return UnitCellLattice.isConsistentWithDesign(self)

    # === MANIPULATION METHODS
    def applyTransF(self, TransF):
        if (isinstance(TransF, ScaleFunc)):
            self._IAD *= TransF.Scale
        UnitCellLattice.applyTransF(self, TransF)

    # === PROPERTY EVALUATION METHODS
    # def isOnLattice(self,P):

    def areNeighbors(self, P1, P2):
        return np.linalg.norm(P2 - P1) <= self.IAD

    def getNeighbors(self, P):
        RefP = self._getConvertToReference(P)
        result = deepcopy(self._RefNeighborsPattern)
        for NeighP in result:
            NeighP += RefP
            self._convertFromReference(NeighP)
        return result

    # === BASIC QUERY METHODS
    @property
    def IAD(self):
        return self._IAD
