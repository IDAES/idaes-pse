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
from math import cos, sin

from ..util.util import areEqual


class TransformFunc(object):
    """ """

    DBL_TOL = 1e-5

    # === PROPERTY EVALUATION METHODS
    @abstractmethod
    def transform(self, P):
        """

        Args:
            P:

        Returns:

        """
        raise NotImplementedError

    @abstractmethod
    def undo(self, P):
        """

        Args:
            P:

        Returns:

        """
        raise NotImplementedError

    def getTransform(self, P):
        """

        Args:
            P:

        Returns:

        """
        result = deepcopy(P)
        self.transform(result)
        return result

    def getUndo(self, P):
        """

        Args:
            P:

        Returns:

        """
        result = deepcopy(P)
        self.undo(result)
        return result

    def __add__(self, other):
        result = CompoundTransformFunc()
        result += self
        result += other
        return result


class ShiftFunc(TransformFunc):
    """ """

    # === STANDARD CONSTRUCTOR
    def __init__(self, Shift):
        self._Shift = Shift

    # === PROPERTY EVALUATION METHODS
    def transform(self, P):
        """

        Args:
            P:

        Returns:

        """
        P += self._Shift

    def undo(self, P):
        """

        Args:
            P:

        Returns:

        """
        P -= self._Shift

    # === BASIC QUERY METHODS
    @property
    def Shift(self):
        """ """
        return self._Shift


class ScaleFunc(TransformFunc):
    """ """

    # === STANDARD CONSTRUCTOR
    def __init__(self, Scale, OriginOfScale=None):
        if isinstance(Scale, np.ndarray):
            self._Scale = Scale
        elif isinstance(Scale, float) or isinstance(Scale, int):
            self._Scale = np.array([Scale, Scale, Scale])
        self._OriginOfScale = OriginOfScale

    # === PROPERTY EVALUATION METHODS
    def transform(self, P):
        """

        Args:
            P:

        Returns:

        """
        if self.OriginOfScale is not None:
            P -= self.OriginOfScale
        P *= self.Scale  # element-wise multiplication
        if self.OriginOfScale is not None:
            P += self.OriginOfScale

    def undo(self, P):
        """

        Args:
            P:

        Returns:

        """
        if self.OriginOfScale is not None:
            P -= self.OriginOfScale
        P /= self.Scale  # element-wise division
        if self.OriginOfScale is not None:
            P += self.OriginOfScale

    # === PROPERTY EVALUATION METHODS
    @property
    def Scale(self):
        """ """
        return self._Scale

    @property
    def OriginOfScale(self):
        """ """
        return self._OriginOfScale

    @property
    def isIsometric(self):
        """ """
        return (self.Scale == self.Scale[0]).all()


class RotateFunc(TransformFunc):
    """ """

    # === STANDARD CONSTRUCTOR
    def __init__(self, RotMat, OriginOfRotation=None):
        self._RotMat = RotMat
        self._RevRotMat = RotMat.transpose()
        self._OriginOfRotation = OriginOfRotation

    # === CONSTRUCTOR - From extrinsic x-y-z series of rotations
    @classmethod
    def fromXYZAngles(cls, ThetaX, ThetaY, ThetaZ, OriginOfRotation=None):
        """

        Args:
            ThetaX: param ThetaY:
            ThetaZ: param OriginOfRotation:  (Default value = None)
            ThetaY:
            OriginOfRotation:  (Default value = None)

        Returns:

        """
        RotMatX = np.array(
            [[1, 0, 0], [0, cos(ThetaX), -sin(ThetaX)], [0, sin(ThetaX), cos(ThetaX)]]
        )
        RotMatY = np.array(
            [[cos(ThetaY), 0, sin(ThetaY)], [0, 1, 0], [-sin(ThetaY), 0, cos(ThetaY)]]
        )
        RotMatZ = np.array(
            [[cos(ThetaZ), -sin(ThetaZ), 0], [sin(ThetaZ), cos(ThetaZ), 0], [0, 0, 1]]
        )
        # import code; code.interact(local=dict(locals(),**globals()));
        return cls(
            np.dot(np.dot(RotMatZ, RotMatY), RotMatX), OriginOfRotation=OriginOfRotation
        )

    @classmethod
    def fromAxisAngle(cls, Axis, Angle, OriginOfRotation=None):
        """

        Args:
            Axis: param Angle:
            OriginOfRotation: Default value = None)
            Angle:

        Returns:

        """
        RotMat = np.zeros((3, 3), dtype=float)
        Axis /= np.linalg.norm(Axis)
        C = cos(Angle)
        S = sin(Angle)
        RotMat[0][0] = C + pow(Axis[0], 2) * (1 - C)
        RotMat[0][1] = Axis[0] * Axis[1] * (1 - C) - Axis[2] * S
        RotMat[0][2] = Axis[0] * Axis[2] * (1 - C) + Axis[1] * S
        RotMat[1][0] = Axis[1] * Axis[0] * (1 - C) + Axis[2] * S
        RotMat[1][1] = C + pow(Axis[1], 2) * (1 - C)
        RotMat[1][2] = Axis[1] * Axis[2] * (1 - C) - Axis[0] * S
        RotMat[2][0] = Axis[2] * Axis[0] * (1 - C) - Axis[1] * S
        RotMat[2][1] = Axis[2] * Axis[1] * (1 - C) + Axis[0] * S
        RotMat[2][2] = C + pow(Axis[2], 2) * (1 - C)
        return cls(RotMat, OriginOfRotation=OriginOfRotation)

    # === PROPERTY EVALUATION METHODS
    def transform(self, P):
        """

        Args:
            P:

        Returns:

        """
        if isinstance(P, np.ndarray) and P.shape == (3, 3):
            # Case of Alignment matrix
            np.dot(self.RotMat, P, out=P)
        else:
            # Case of point (np.array)
            if self.OriginOfRotation is not None:
                P -= self.OriginOfRotation
            np.dot(self.RotMat, P, out=P)  # matrix multiplication
            if self.OriginOfRotation is not None:
                P += self.OriginOfRotation

    def transformDirection(self, V):
        """

        Args:
            V:

        Returns:

        """
        np.dot(self.RotMat, V, out=V)

    def undo(self, P):
        """

        Args:
            P:

        Returns:

        """
        if self.OriginOfRotation is not None:
            P -= self.OriginOfRotation
        np.dot(self.RevRotMat, P, out=P)
        if self.OriginOfRotation is not None:
            P += self.OriginOfRotation

    def undoDirection(self, V):
        """

        Args:
            V:

        Returns:

        """
        np.dot(self.RevRotMat, V, out=V)

    # === BASIC QUERY METHODS
    @property
    def RotMat(self):
        """ """
        return self._RotMat

    @property
    def RevRotMat(self):
        """ """
        return self._RevRotMat

    @property
    def OriginOfRotation(self):
        """ """
        return self._OriginOfRotation


class ReflectFunc(TransformFunc):
    """ """

    # === STANDARD CONSTRUCTOR
    def __init__(self, PlaneNorm, PlaneRefPoint):
        assert areEqual(np.linalg.norm(PlaneNorm), 1.0, TransformFunc.DBL_TOL)
        self._PlaneNorm = PlaneNorm
        self._PlaneRefPoint = PlaneRefPoint

    # === CONSTRUCTOR - Reflection across X axis
    @classmethod
    def fromPoints(cls, P0, P1, P2):
        """

        Args:
            P0: param P1:
            P2:
            P1:

        Returns:

        """
        TmpPlaneNorm = np.cross(P1 - P0, P2 - P0)
        TmpPlaneNorm /= np.linalg.norm(TmpPlaneNorm)
        return cls(TmpPlaneNorm, P0)

    @classmethod
    def acrossX(cls):
        """ """
        return cls(np.array([1, 0, 0], dtype=float), np.array([0, 0, 0], dtype=float))

    @classmethod
    def acrossY(cls):
        """ """
        return cls(np.array([0, 1, 0], dtype=float), np.array([0, 0, 0], dtype=float))

    @classmethod
    def acrossZ(cls):
        """ """
        return cls(np.array([0, 0, 1], dtype=float), np.array([0, 0, 0], dtype=float))

    # === PROPERTY EVALUATION METHODS
    def transform(self, P):
        """

        Args:
            P:

        Returns:

        """
        if type(P) is np.ndarray and P.shape == (3, 3):
            # Case of Alignment matrix
            for AlignmentAxis in P:  # i.e., apply to each vector, but do not shift
                Vperp = np.inner(AlignmentAxis, self.PlaneNorm) * self.PlaneNorm
                AlignmentAxis -= 2 * Vperp
        else:
            # Case of Point
            V = P - self.PlaneRefPoint
            Vperp = np.inner(V, self.PlaneNorm) * self.PlaneNorm
            P -= 2 * Vperp

    def undo(self, P):
        """

        Args:
            P:

        Returns:

        """
        self.transform(P)  # same as transform

    # === BASIC QUERY METHODS
    @property
    def PlaneNorm(self):
        """ """
        return self._PlaneNorm

    @property
    def PlaneRefPoint(self):
        """ """
        return self._PlaneRefPoint


class CompoundTransformFunc(TransformFunc):
    """ """

    # === STANDARD CONSTRUCTOR
    def __init__(self, args=None):
        self._TransFs = args if args is not None else []

    # === PROPERTY EVALUATION METHODS
    def transform(self, P):
        """

        Args:
            P:

        Returns:

        """
        for TransF in self.TransFs:
            TransF.transform(P)

    def undo(self, P):
        """

        Args:
            P:

        Returns:

        """
        for TransF in reversed(self.TransFs):
            TransF.undo(P)

    # === BASIC QUERY METHODS
    @property
    def TransFs(self):
        """ """
        return self._TransFs

    def __iadd__(self, other):
        if isinstance(other, CompoundTransformFunc):
            self._TransFs.extend(other.TransFs)
        elif isinstance(other, TransformFunc):
            self._TransFs.append(other)
        else:
            raise ValueError("Decide how to handle this case!")
        return self
