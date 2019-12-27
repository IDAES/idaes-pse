from abc import abstractmethod
from copy import deepcopy
import numpy as np

from ..transform_func import TransformFunc, ShiftFunc, ScaleFunc, RotateFunc, ReflectFunc


class Lattice(object):
    '''
    # === AUXILIARY METHODS
    @abstractmethod
    def _isValidScanPoint(self,P):
        raise NotImplementedError

    def _setupScan(self,ScanMin,ScanMax):
        self._ScanMin = ScanMin.astype(int)+np.array([-1,-1,-1],dtype=int)
        self._ScanMax = ScanMax.astype(int)+np.array([ 1, 1, 1],dtype=int)
        self._ScanP = deepcopy(self._ScanMin)
        self._ScanP[0] -= 1
        #print('Setting up ScanMin={}'.format(self._ScanMin))
        #print('Setting up ScanMax={}'.format(self._ScanMax))
        #print('Setting up ScanP={}'.format(self._ScanP))

    def _getNextPoint(self):
        #import code; code.interact(local=dict(locals(),**globals()));
        while(self._ScanP[2]<=self._ScanMax[2]):
            while(self._ScanP[1]<=self._ScanMax[1]):
                while(self._ScanP[0]<=self._ScanMax[0]):
                    self._ScanP[0] += 1
                    if(self._ScanP[0]<=self._ScanMax[0] and
                       self._ScanP[1]<=self._ScanMax[1] and
                       self._ScanP[2]<=self._ScanMax[2] and 
                       self._isValidScanPoint(self._ScanP)):
                        #print('Returning P={}'.format(self._getConvertFromReference(self._ScanP.astype(float))))
                        return self._getConvertFromReference(self._ScanP.astype(float))
                self._ScanP[0] = self._ScanMin[0]
                self._ScanP[1] += 1
            self._ScanP[1] = self._ScanMin[1]
            self._ScanP[2] += 1
        self._ScanP[2] = self._ScanMin[2]
        return None # only invoked when ScanP is reset to ScanMin
    '''

    # === DEFAULT CONSTRUCTOR
    def __init__(self):
        self._TransformFuncs = []

    # === MANIPULATION METHODS
    def applyTransF(self, TransF):
        if (isinstance(TransF, TransformFunc)):
            self._TransformFuncs.append(TransF)
        else:
            raise TypeError

    def shift(Shift):
        if (type(Shift) is ShiftFunc):
            self.applyTransF(Shift)
        elif (type(Shift) is np.ndarray):
            self.applyTransF(ShiftFunc(Shift))
        else:
            raise TypeError

    def scale(Scale, OriginOfScale=None):
        if (type(Scale) is ScaleFunc):
            self.applyTransF(Scale)
        elif (type(Scale) is np.ndarray):
            self.applyTransF(ScaleFunc(Scale, OriginOfScale))
        else:
            raise TypeError

    def rotate(Rotation, OriginOfRotation=None):
        if (type(Rotation) is RotateFunc):
            self.applyTransF(Rotation)
        elif (type(Rotation) is np.ndarray):
            self.applyTransF(RotateFunc(Rotation, OriginOfRotation))
        else:
            raise TypeError

    def reflect(Reflection):
        if (type(Reflection) is ReflectFunc):
            self.applyTransF(Reflection)
        else:
            raise TypeError

    # === PROPERTY EVALUATION METHODS
    @abstractmethod
    def isOnLattice(self, P):
        raise NotImplementedError

    @abstractmethod
    def areNeighbors(self, P1, P2):
        raise NotImplementedError

    @abstractmethod
    def getNeighbors(self, P):
        raise NotImplementedError

    def _convertFromReference(self, P):
        for TransF in self._TransformFuncs:
            TransF.transform(P)

    def _convertToReference(self, P):
        for TransF in reversed(self._TransformFuncs):
            TransF.undo(P)

    def _getConvertFromReference(self, P):
        result = deepcopy(P)
        self._convertFromReference(result)
        return result

    def _getConvertToReference(self, P):
        result = deepcopy(P)
        self._convertToReference(result)
        return result
