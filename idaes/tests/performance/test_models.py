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

from idaes.core.solvers import get_solver

from idaes.models.properties.modular_properties.examples.tests.test_HC_PR import (
    HC_PR_Model,
)
from idaes.models.unit_models.tests.test_heat_exchanger_1D import HX1D_Model


# Get default solver for testing
solver = get_solver()


class TestModel(unittest.TestCase):
    """
    This class is intended to be a generic template for common performance testing activities.
    It is not intended to be universal however, and test developers should feel free to create their
    own test classes if required (e.g. to test performance of different aspects of model performance
    or to get more granular information).
    """

    def recordData(self, name, value):
        """A method for recording data associated with a test.

        This method is only meaningful when running this TestCase with
        'nose', using the TestData plugin.

        """
        tmp = getattr(self, "testdata", None)
        if not tmp is None:
            tmp[name] = value

    def _run_test(self, test_case, consistent_units=True):
        """
        General method for executing common performance test runs.

        This method takes a test case class which is expected to implement three methods;
        i) build_model, ii) initialize_model and iii) solve_model. Execution times is measured for
        each of these methods, along with time taken to check unit consistency.

        Args:
            test_case - class containing methods required to construct, initialize and solve test case.
            consistent_units - bool indicating whether unit consistency check should return True (default) or False.
        """
        # Build model and record execution time
        gc.collect()
        timer = TicTocTimer()
        model = test_case.build_model()
        self.recordData("build model", timer.toc("build model"))

        # Check unit consistency and record execution time
        gc.collect()
        timer.tic(None)
        try:
            assert_units_consistent(model)

            if not consistent_units:
                # If consistent_units is False we expected this to fail, so raise Exception
                raise Exception("Unit consistency test did not fail.")
        except:
            # If consistent_units is False, we can ignore this exception
            if not consistent_units:
                pass
            else:
                raise
        finally:
            self.recordData("unit consistency", timer.toc("unit consistency"))

        # Initialize model and record execution time
        gc.collect()
        timer.tic(None)
        test_case.initialize_model(model)
        self.recordData("initialize", timer.toc("initialize"))

        # Solve model and record execution time
        gc.collect()
        timer.tic(None)
        try:
            test_case.solve_model(model)
        except AttributeError:
            solver.solve(model)
        self.recordData("final solve", timer.toc("final solve"))


@pytest.mark.performance
class TestIdaesPerformance(TestModel):
    def test_hydrocarbon_PR_properties(self):
        self._run_test(HC_PR_Model)

    def test_heat_exchanger_1D_IAPWS(self):
        self._run_test(HX1D_Model)
