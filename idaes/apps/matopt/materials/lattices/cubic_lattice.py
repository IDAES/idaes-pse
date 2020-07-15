##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
import numpy as np
from copy import deepcopy

from ..geometry import Parallelepiped
from ..transform_func import ScaleFunc
from .unit_cell_lattice import UnitCell, UnitCellLattice
from ..tiling import CubicTiling


class CubicLattice(UnitCellLattice):
    RefIAD = 1

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

    # === MANIPULATION METHODS
    def applyTransF(self, TransF):
        if isinstance(TransF, ScaleFunc):
            self._IAD *= TransF.Scale
        UnitCellLattice.applyTransF(self, TransF)

    # === PROPERTY EVALUATION METHODS
    # NOTE: This method is inherited from UnitCellLattice
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
