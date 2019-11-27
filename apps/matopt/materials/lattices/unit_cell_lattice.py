
from abc import abstractmethod
from copy import deepcopy
import numpy as np

from ...util.util import myArrayEq
from .lattice import Lattice
from ..transform_func import TransformFunc

class UnitCell(object):
    DBL_TOL = 1e-5

    # === STANDARD CONSTRUCTOR
    def __init__(self,argTiling,FracPositions):
        self._Tiling = argTiling
        self._FracPositions = FracPositions

    # === ASSERTION OF CLASS DESIGN
    def isConsistentWithDesign(self):
        for Position in self.FracPositions:
            for coord in Position:
                if(coord < 0.0 - UnitCell.DBL_TOL or
                   coord > 1.0 + UnitCell.DBL_TOL):
                    return False
        return self.Tiling.isConsistentWithDesign()

    # === MANIPULATION METHODS
    def addPosition(self,P):
        self._FracPositions.append(P)

    def applyTransF(self,TransF):
        if(isinstance(TransF,TransformFunc)):
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

    def convertToFrac(self,P,blnDiscardIntPart=True):
        P[:] = self.Tiling.getFractionalCoords(P,
                                               blnRelativeToCenter=False,
                                               blnRoundInside=True,
                                               blnPreferZero=True)
        if(blnDiscardIntPart):
            P -= P.astype(int)
            '''
            P = np.array([P[0]%1,
                          P[1]%1,
                          P[2]%1],dtype=float)
            '''
    def getConvertToFrac(self,P,blnDiscardIntPart=True):
        result = deepcopy(P)
        self.convertToFrac(result,blnDiscardIntPart)
        return result

    def recenter(self,P):
        while(np.inner(P-PositiveXYZCorner, Nx) > UnitCell.DBL_TOL): P -= Vx
        while(np.inner(P-NegativeXYZCorner,-Nx) > UnitCell.DBL_TOL): P += Vx
        while(np.inner(P-PositiveXYZCorner, Ny) > UnitCell.DBL_TOL): P -= Vy
        while(np.inner(P-NegativeXYZCorner,-Ny) > UnitCell.DBL_TOL): P += Vy
        while(np.inner(P-PositiveXYZCorner, Nz) > UnitCell.DBL_TOL): P -= Vz
        while(np.inner(P-NegativeXYZCorner,-Nz) > UnitCell.DBL_TOL): P += Vz

    def getRecentered(self,P):
        result = deepcopy(P)
        self.recenter(result)
        return result

    def getPointType(self,P):
        PFrac = self.getConvertToFrac(P,blnDiscardIntPart=True)
        for i,FracPosition in enumerate(self.FracPositions):
            if(myArrayEq(FracPosition,PFrac,UnitCell.DBL_TOL)):
                return i
        return None


class UnitCellLattice(Lattice):
    # === AUXILIARY METHODS
    def ScanRef(self,RefScanMin,RefScanMax):
        #print('{} {}'.format(RefScanMin,RefScanMax))
        for n1 in range(RefScanMin[0],RefScanMax[0]+1):
            for n2 in range(RefScanMin[1],RefScanMax[1]+1):
                for n3 in range(RefScanMin[2],RefScanMax[2]+1):
                    for FracPart in self.RefUnitCell.FracPositions:
                        yield (n1*self.RefUnitCell.Tiling.TileShape.Vx+
                               n2*self.RefUnitCell.Tiling.TileShape.Vy+
                               n3*self.RefUnitCell.Tiling.TileShape.Vz+
                               FracPart)

    def Scan(self,argPolyhedron):
        RefScanMin = self._getConvertToReference(argPolyhedron.V[0])
        RefScanMax = self._getConvertToReference(argPolyhedron.V[0])
        for v in argPolyhedron.V:
            RefScanMin = np.minimum(RefScanMin,self._getConvertToReference(v))
            RefScanMax = np.maximum(RefScanMax,self._getConvertToReference(v))
        # Add some padding to these numbers to handle integer rounding 
        #import code; code.interact(local=dict(locals(),**globals()));
        RefScanMin = RefScanMin.astype(int)+np.array([-1,-1,-1],dtype=int)
        RefScanMax = RefScanMax.astype(int)+np.array([ 1, 1, 1],dtype=int) 
        for RefP in self.ScanRef(RefScanMin,RefScanMax):
            yield self._getConvertFromReference(RefP)


    '''
    def _getNextRefPoint(self):
        #import code; code.interact(local=dict(locals(),**globals()));
        while(self._InterCellScanP[2]<=self._ScanMax[2]):
            while(self._InterCellScanP[1]<=self._ScanMax[1]):
                while(self._InterCellScanP[0]<=self._ScanMax[0]):
                    if(self._InterCellScanP[0]<=self._ScanMax[0] and
                       self._InterCellScanP[1]<=self._ScanMax[1] and
                       self._InterCellScanP[2]<=self._ScanMax[2]):
                        while(self._IntraCellScanI < len(self.RefUnitCell.FracPositions)):
                            result = (self._InterCellScanP+
                                      self.RefUnitCell.FracPositions[self._IntraCellScanI])
                            self._IntraCellScanI += 1
                            return result
                        self._InterCellScanP[0] += 1
                        self._IntraCellScanI = 0
                self._InterCellScanP[0] = self._ScanMin[0]
                self._InterCellScanP[1] += 1
            self._InterCellScanP[1] = self._ScanMin[1]
            self._InterCellScanP[2] += 1
        self._InterCellScanP[2] = self._ScanMin[2]
        return None # only invoked when ScanP is reset to ScanMin
    '''

    # === STANDARD CONSTRUCTOR
    def __init__(self,RefUnitCell):
        self._RefUnitCell = RefUnitCell
        Lattice.__init__(self)
        assert(self.isConsistentWithDesign())

    # === ASSERTION OF CLASS DESIGN
    def isConsistentWithDesign(self):
        # TODO: Add some real checks here
        return True

    # === MANIPULATION METHODS

    # === PROPERTY EVALUATION METHODS
    def isOnLattice(self,P):
        return self.RefUnitCell.getPointType(Lattice._getConvertToReference(self,P)) is not None

    @abstractmethod
    def areNeighbors(self,P1,P2):
        raise NotImplementedError

    @abstractmethod
    def getNeighbors(self,P):
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
