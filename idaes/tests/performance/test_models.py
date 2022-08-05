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
"""
Author: Andrew Lee
"""
import pytest
import gc

import pyomo.common.unittest as unittest
from pyomo.util.check_units import assert_units_consistent
from pyomo.common.timing import TicTocTimer

from idaes.models.properties.modular_properties.examples.tests import test_HC_PR


class TestModel(unittest.TestCase):
    def recordData(self, name, value):
        """A method for recording data associated with a test.

        This method is only meaningful when running this TestCase with
        'nose', using the TestData plugin.

        """
        tmp = getattr(self, "testdata", None)
        if not tmp is None:
            tmp[name] = value

    def _run_test(self, model_lib, init_model, solve_model):
        gc.collect()
        timer = TicTocTimer()

        model = model_lib()
        self.recordData("build model", timer.toc("build model"))

        gc.collect()
        timer.tic(None)
        assert_units_consistent(model)
        self.recordData("unit consistency", timer.toc("unit consistency"))

        gc.collect()
        timer.tic(None)
        init_model(model)
        self.recordData("initialize", timer.toc("initialize"))

        gc.collect()
        timer.tic(None)
        solve_model(model)
        self.recordData("final solve", timer.toc("final solve"))


@pytest.mark.performance
class TestMisc(TestModel):
    def test_hydrocarbon_PR_properties(self):
        self._run_test(
            test_HC_PR.build_model, test_HC_PR.initialize_model, test_HC_PR.solve_model
        )
