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

from ..transform_func import (
    TransformFunc,
    ShiftFunc,
    ScaleFunc,
    RotateFunc,
    ReflectFunc,
)


class Lattice(object):
    """A class used to represent crystal lattice locations.

    The class encodes methods for determining which Cartesian coordinates to
    consider as sites on an infinite crystal lattice. A ``Lattice`` can be constructed from
    a point on the lattice (i.e., a shift from the origin), an alignment (i.e., rotation from a
    nominal orientation), and appropriate scaling factors. With these attributes, we generally
    support the translation, rotation, and rescaling of lattices. Additionally, ``Lattice`` objects
    include a method for determining which sites should be considered neighbors.
    """

    # === DEFAULT CONSTRUCTOR
    def __init__(self):
        self._TransformFuncs = []

    # === MANIPULATION METHODS
    def applyTransF(self, TransF):
        if isinstance(TransF, TransformFunc):
            self._TransformFuncs.append(TransF)
        else:
            raise TypeError

    def shift(self, Shift):
        if type(Shift) is ShiftFunc:
            self.applyTransF(Shift)
        elif type(Shift) is np.ndarray:
            self.applyTransF(ShiftFunc(Shift))
        else:
            raise TypeError

    def scale(self, Scale, OriginOfScale=None):
        if type(Scale) is ScaleFunc:
            self.applyTransF(Scale)
        elif type(Scale) is np.ndarray:
            self.applyTransF(ScaleFunc(Scale, OriginOfScale))
        else:
            raise TypeError

    def rotate(self, Rotation, OriginOfRotation=None):
        if type(Rotation) is RotateFunc:
            self.applyTransF(Rotation)
        elif type(Rotation) is np.ndarray:
            self.applyTransF(RotateFunc(Rotation, OriginOfRotation))
        else:
            raise TypeError

    def reflect(self, Reflection):
        if type(Reflection) is ReflectFunc:
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
    def getNeighbors(self, P, layer):
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
