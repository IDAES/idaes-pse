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
import gc

import pyomo.common.unittest as unittest
from pyomo.util.check_units import assert_units_consistent
from pyomo.common.timing import TicTocTimer
from pyomo.environ import assert_optimal_termination

from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom

# Get default solver for testing
solver = get_solver()


class PerformanceBaseClass:  # unittest.TestCase):
    """
    Base class for running the IDAES Performance Test Suite.

    This class is intended to be used as a basis for constructing tests to be run as part of the IDAES
    performance testing suite, and should implement the following methods:

    * build_model - should return a complete instance of the model being tested.
    * test_performance - a test method for profiling the model of interest which logs some performance data.
        A standard method is provided which logs some common performance data, but developers should overload this
        as required.
    * initialize_model - (optional) a method to initialize the model of interest. This method is called as part
        of the standard test_performance method.
    * solve_model - (optional) a method to solve the model of interest. This method is called as part
        of the standard test_performance method, and a standard method is provided by the base class.

    Developers may add additional test methods to the derived class to run additional performance tests if desired.
    These methods should ass start with 'test_' in order ot be discovered by pytest.

    Usage:

    Performance tests should inherit from this based class and decorate the derived class with the 'performance'
    pytest mark. This will automatically categorize all tests in this class as performance tests which will only be run
    when explicitly called for.
    """

    def recordData(self, name, value):
        """
        A method for recording data associated with a test.
        """
        tmp = getattr(self, "testdata", None)
        if not tmp is None:
            tmp[name] = value

    def build_model(self):
        """
        Return a fully constructed version of the model of interest. Developers must overload this method.

        Args:
            None

        Returns:
            None
        """
        raise NotImplementedError(
            "Test class has not implemented a build_model method."
        )

    def initialize_model(self, model):
        """
        Initialize provided model. Developers should overload this method to provide model specific initialization
        instructions.

        Args:
            model - constructed model object to be initialized

        Returns:
            None
        """
        raise NotImplementedError(
            "Test class has not implemented an initialize_model method."
        )

    def solve_model(self, model):
        """
        Solve provided model. Developers should overload this if necessary to customize solver and arguments.

        Args:
            model - constructed model object to be solved

        Returns:
            None
        """
        results = solver.solve(model)

        assert_optimal_termination(results)

    def test_performance(self):
        """
        Standard method for executing performance test runs.

        This method calls the build_model, initialize_model and solve_model methods thus developed must provide these.

        This method logs the following information:

        1. Model construction time
        2. Time required to assert unit consistency
        3. Model initialization time
        4. Model solution time (post-initialization)

        Developers may overload this method to customize the information recorded. Note that dynamic models currently
        fail unit consistency checks, and thus method will fail if used on a dynamic model (for now).

        Args:
            None
        """
        # Build model and record execution time
        gc.collect()
        timer = TicTocTimer()
        model = self.build_model()
        self.recordData("build model", timer.toc("build model"))
        assert degrees_of_freedom(model) == 0

        # Check unit consistency and record execution time
        gc.collect()
        timer.tic(None)
        assert_units_consistent(model)
        self.recordData("unit consistency", timer.toc("unit consistency"))

        # Initialize model and record execution time
        gc.collect()
        timer.tic(None)
        self.initialize_model(model)
        self.recordData("initialize", timer.toc("initialize"))

        # Solve model and record execution time
        gc.collect()
        timer.tic(None)
        self.solve_model(model)
        self.recordData("final solve", timer.toc("final solve"))
