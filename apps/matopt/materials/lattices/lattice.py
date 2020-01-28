from abc import abstractmethod
from copy import deepcopy
import numpy as np

from ..transform_func import TransformFunc, ShiftFunc, ScaleFunc, RotateFunc, ReflectFunc


class Lattice(object):
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
