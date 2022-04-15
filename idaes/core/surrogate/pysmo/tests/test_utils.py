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
#  ___________________________________________________________________________
#
#  Pyomo: Python Optimization Modeling Objects
#  Copyright 2017 National Technology and Engineering Solutions of Sandia, LLC
#  Under the terms of Contract DE-NA0003525 with National Technology and
#  Engineering Solutions of Sandia, LLC, the U.S. Government retains certain
#  rights in this software.
#  This software is distributed under the 3-clause BSD License.
#  ___________________________________________________________________________

import pyomo.common.unittest as unittest
import pytest
from idaes.core.surrogate.pysmo.utils import NumpyEvaluator
from pyomo.environ import ConcreteModel, Param, Var, value, sin, atan, atanh
from pyomo.core import ComponentMap

try:
    import numpy as np

    _numpy_available = True
except ImportError:
    _numpy_available = False


@unittest.skipIf(not _numpy_available, "Test requires numpy")
class TestNumpyEvaluator:
    @pytest.mark.unit
    def test_eval_numpy(self):
        m = ConcreteModel()
        m.p = Param([1, 2], mutable=True)
        m.x = Var()

        data = np.array([[0, -1, 2], [0.1, 0.2, 0.3], [4, 5, 6]])
        cMap = ComponentMap()
        cMap[m.p[1]] = data[0]
        cMap[m.p[2]] = data[1]
        cMap[m.x] = data[2]

        npe = NumpyEvaluator(cMap)

        result = npe.walk_expression(sin(m.x))
        assert pytest.approx(result[0], rel=1e-12) == sin(4)
        assert pytest.approx(result[1], rel=1e-12) == sin(5)
        assert pytest.approx(result[2], rel=1e-12) == sin(6)

        result = npe.walk_expression(abs(m.x * m.p[1] - m.p[2]))
        assert pytest.approx(result[0], rel=1e-12) == 0.1
        assert pytest.approx(result[1], rel=1e-12) == -((-1 * 5) - 0.2)
        assert pytest.approx(result[2], rel=1e-12) == (2 * 6 - 0.3)

        result = npe.walk_expression(atan(m.x))
        assert pytest.approx(result[0], rel=1e-12) == atan(4)
        assert pytest.approx(result[1], rel=1e-12) == atan(5)
        assert pytest.approx(result[2], rel=1e-12) == atan(6)

        result = npe.walk_expression(atanh(m.p[2]))
        assert pytest.approx(result[0], rel=1e-12) == atanh(0.1)
        assert pytest.approx(result[1], rel=1e-12) == atanh(0.2)
        assert pytest.approx(result[2], rel=1e-12) == atanh(0.3)

    @pytest.mark.unit
    def test_eval_constant(self):
        m = ConcreteModel()
        m.p = Param([1, 2], mutable=True)
        m.x = Var(initialize=0.25)

        cMap = ComponentMap()
        cMap[m.p[1]] = 2
        cMap[m.p[2]] = 4

        npe = NumpyEvaluator(cMap)

        expr = m.p[1] + m.p[2] + m.x + 0.5
        assert npe.walk_expression(expr) == pytest.approx(6.75, rel=1e-12)

        m.p[1] = 2
        m.p[2] = 4
        assert value(expr) == pytest.approx(6.75, rel=1e-12)


if __name__ == "__main__":
    pytest.main()
