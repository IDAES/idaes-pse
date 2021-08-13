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
Tests for Alampy SurrogateModelTrainer
"""
import pytest
import numpy as np
import io

from idaes.surrogate.alamopy_new import Alamopy, Modelers, Screener


class TestAlamopyWriter():
    @pytest.fixture
    def alm_obj(self):
        alm_obj = Alamopy()

        alm_obj._n_inputs = 2
        alm_obj._n_outputs = 1

        alm_obj._input_labels = ["x1", "x2"]
        alm_obj._output_labels = ["z1"]

        alm_obj._input_min = [0, 0]
        alm_obj._input_max = [5, 10]

        alm_obj._rdata_in = np.array([[1, 2, 3, 4], [5, 6, 7, 8]])
        alm_obj._rdata_out = np.array([10, 20, 30, 40], ndmin=2)

        return alm_obj

    @pytest.mark.unit
    def test_default(self, alm_obj):
        stream = io.StringIO()
        alm_obj._write_alamo_input_to_stream(stream=stream)

        assert stream.getvalue() == (
            "# IDAES Alamopy input file\n"
            "NINPUTS 2\n"
            "NOUTPUTS 1\n"
            "XLABELS x1 x2\n"
            "ZLABELS z1\n"
            "XMIN 0 0\n"
            "XMAX 5 10\n"
            "NDATA 4\n"
            "NVALDATA 0\n\n"
            "linfcns 1\n"
            "constant 1\n"
            "maxtime 1000.0\n"
            "numlimitbasis 1\n\n"
            "TRACE 1\n\n"
            "BEGIN_DATA\n"
            "1 5 10\n"
            "2 6 20\n"
            "3 7 30\n"
            "4 8 40\n"
            "END_DATA\n")

    @pytest.mark.unit
    def test_trace(self, alm_obj):
        stream = io.StringIO()
        alm_obj._write_alamo_input_to_stream(
            stream=stream, trace_fname="foo.bar")

        assert stream.getvalue() == (
            "# IDAES Alamopy input file\n"
            "NINPUTS 2\n"
            "NOUTPUTS 1\n"
            "XLABELS x1 x2\n"
            "ZLABELS z1\n"
            "XMIN 0 0\n"
            "XMAX 5 10\n"
            "NDATA 4\n"
            "NVALDATA 0\n\n"
            "linfcns 1\n"
            "constant 1\n"
            "maxtime 1000.0\n"
            "numlimitbasis 1\n\n"
            "TRACE 1\n"
            "TRACEFNAME foo.bar\n"
            "BEGIN_DATA\n"
            "1 5 10\n"
            "2 6 20\n"
            "3 7 30\n"
            "4 8 40\n"
            "END_DATA\n")

    @pytest.mark.unit
    def test_vdata(self, alm_obj):
        alm_obj._vdata_in = np.array([[2.5], [6.5]])
        alm_obj._vdata_out = np.array([25], ndmin=2)

        stream = io.StringIO()
        alm_obj._write_alamo_input_to_stream(stream=stream)

        assert stream.getvalue() == (
            "# IDAES Alamopy input file\n"
            "NINPUTS 2\n"
            "NOUTPUTS 1\n"
            "XLABELS x1 x2\n"
            "ZLABELS z1\n"
            "XMIN 0 0\n"
            "XMAX 5 10\n"
            "NDATA 4\n"
            "NVALDATA 1\n\n"
            "linfcns 1\n"
            "constant 1\n"
            "maxtime 1000.0\n"
            "numlimitbasis 1\n\n"
            "TRACE 1\n\n"
            "BEGIN_DATA\n"
            "1 5 10\n"
            "2 6 20\n"
            "3 7 30\n"
            "4 8 40\n"
            "END_DATA\n"
            "\nBEGIN_VALDATA\n"
            "2.5 6.5 25\n"
            "END_VALDATA\n")

    @pytest.mark.unit
    def test_full_config(self, alm_obj):
        alm_obj.config.xfactor = [1.2, 1.6]
        alm_obj.config.xscaling = True
        alm_obj.config.scalez = True

        alm_obj.config.monomialpower = [1, 2]
        alm_obj.config.multi2power = [1, 2, 3]
        alm_obj.config.multi3power = [1, 2, 3, 4]
        alm_obj.config.ratiopower = [1, 2, 3, 4, 5]

        alm_obj.config.constant = True
        alm_obj.config.linfcns = False
        alm_obj.config.expfcns = True
        alm_obj.config.logfcns = False
        alm_obj.config.sinfcns = True
        alm_obj.config.cosfcns = False
        alm_obj.config.grbfcns = True
        alm_obj.config.rbfparam = 7

        alm_obj.config.modeler = Modelers.AICc
        alm_obj.config.builder = True
        alm_obj.config.backstepper = False

        stream = io.StringIO()
        alm_obj._write_alamo_input_to_stream(stream=stream)

        assert stream.getvalue() == (
            "# IDAES Alamopy input file\n"
            "NINPUTS 2\n"
            "NOUTPUTS 1\n"
            "XLABELS x1 x2\n"
            "ZLABELS z1\n"
            "XMIN 0 0\n"
            "XMAX 5 10\n"
            "NDATA 4\n"
            "NVALDATA 0\n\n"
            "MONO 2\n"
            "MULTI2 3\n"
            "MULTI3 4\n"
            "RATIOS 5\n"
            "xfactor 1.2 1.6\n"
            "xscaling 1\n"
            "scalez 1\n"
            "monomialpower 1 2\n"
            "multi2power 1 2 3\n"
            "multi3power 1 2 3 4\n"
            "ratiopower 1 2 3 4 5\n"
            "expfcns 1\n"
            "linfcns 0\n"
            "logfcns 0\n"
            "sinfcns 1\n"
            "cosfcns 0\n"
            "constant 1\n"
            "grbfcns 1\n"
            "rbfparam 7.0\n"
            "modeler 3\n"
            "builder 1\n"
            "backstepper 0\n"
            "maxtime 1000.0\n"
            "numlimitbasis 1\n\n"
            "TRACE 1\n\n"
            "BEGIN_DATA\n"
            "1 5 10\n"
            "2 6 20\n"
            "3 7 30\n"
            "4 8 40\n"
            "END_DATA\n")
        assert False

    # CONFIG.declare('convpen', ConfigValue(
    #     default=None,
    #     domain=float,
    #     description="Convex penalty term to use if Modeler == SSEP or MADp.",
    #     doc="When MODELER is set to 6 or 8, a penalty consisting of the sum "
    #     "of square errors (SSEP) or the maximum absolute error (MADp) and a "
    #     "term penalizing model size is used for model building. The size of "
    #     "the model is weighted by convpen. If convpen=0, this metric reduces "
    #     "to the classical sum of square errors (SSEP) or the maximum absolute "
    #     "deviation (MADp)."))
    # CONFIG.declare('screener', ConfigValue(
    #     default=None,
    #     domain=In(Screener),
    #     description="Regularization method used to reduce the number of "
    #     "potential basis functions before optimization. Must be instance of "
    #     "Screener Enum."))
    # CONFIG.declare('ncvf', ConfigValue(
    #     default=None,
    #     domain=int,
    #     description="Number of folds to be used for cross validation by the "
    #     "lasso screener. ALAMO will use a two-fold validation if fewer than "
    #     "10 data points are available. NCVF must be a nonnegative integer."))
    # CONFIG.declare('sismult', ConfigValue(
    #     default=None,
    #     domain=int,
    #     description="This parameter must be non-negative and is used to "
    #     "determine the number of basis functions retained by the SIS "
    #     "screener. The number of basis functions retained equals the floor "
    #     "of SSISmult n ln(n), where n is the number of measurements "
    #     "available at the current ALAMO iteration."))

    # CONFIG.declare('maxiter', ConfigValue(
    #     default=None,
    #     domain=int,
    #     description="Maximum number of ALAMO iterations. 1 = no adaptive "
    #     "sampling, 0 = no limit.",
    #     doc="Maximum number of ALAMO iterations. Each iteration begins with "
    #     "a model-building step. An adaptive sampling step follows if maxiter "
    #     "does not equal 1. If maxiter is set to a number less than or equal "
    #     "to 0, ALAMO will enforce no limit on the number of iterations."))
    # CONFIG.declare('maxtime', ConfigValue(
    #     default=1000,
    #     domain=float,
    #     description="Maximum total execution time allowed in seconds. "
    #     "Default = 1000.",
    #     doc="Maximum total execution time allowed in seconds. This time "
    #     "includes all steps of the algorithm, including time to read problem, "
    #     "preprocess data, solve optimization subproblems, and print results."))
    # CONFIG.declare('datalimitterms', ConfigValue(
    #     default=None,
    #     domain=In([True, False]),
    #     description="Limit model terms to number of measurements.",
    #     doc="If True, ALAMO will limit the number of terms in the model to be "
    #     "no more than the number of data measurements; otherwise, no limit "
    #     "based on the number of data measurements will be placed. The user "
    #     "may provide an additional limit on the number of terms in the model "
    #     "through the maxterms and minterms options."))
    # CONFIG.declare('maxterms', ConfigValue(
    #     default=None,
    #     domain=list_of_ints,
    #     description="List of maximum number of model terms to per output.",
    #     doc="Row vector of maximum terms allowed in the modeling of output "
    #     "variables. One per output variable, space separated. A −1 signals "
    #     "that no limit is imposed."))
    # CONFIG.declare('minterms', ConfigValue(
    #     default=None,
    #     domain=list_of_ints,
    #     description="List of minimum number of model terms to per output.",
    #     doc="Row vector of minimum terms required in the modeling of output "
    #     "variables. One per output variable, space separated. A 0 signals "
    #     "that no limit is imposed."))
    # CONFIG.declare('numlimitbasis', ConfigValue(
    #     default=True,  # default this to true to avoid numerical issues
    #     domain=In([True, False]),
    #     description="Eliminate infeasible basis functions. Default = True",
    #     doc="If True, ALAMO will eliminate basis functions that are not "
    #     "numerically acceptable (e.g., log(x) will be eliminated if x may be "
    #     "negative); otherwise, no limit based on the number of data "
    #     "measurements will be placed. The user may provide additional limits "
    #     "on the the type and number of selected basis functions through the "
    #     "options exclude and groupcon."))
    # CONFIG.declare('exclude', ConfigValue(
    #     default=None,
    #     domain=list_of_ints,
    #     description="List of inputs to exclude during building,",
    #     doc="Row vector of 0/1 flags that specify which input variables, if "
    #     "any, ALAMO should exclude during the model building process. All "
    #     "input variables must be present in the data but ALAMO will not "
    #     "include basis functions that involve input variables for which "
    #     "exclude equals 1. This feature does not apply to custom basis "
    #     "functions or RBFs."))
    # CONFIG.declare('ignore', ConfigValue(
    #     default=None,
    #     domain=list_of_ints,
    #     description="List of outputs to ignore during building.",
    #     doc="Row vector of 0/1 flags that specify which output variables, "
    #     "if any, ALAMO should ignore. All output variables must be present in "
    #     "the data but ALAMO does not model output variables for which ignore "
    #     "equals 1."))
    # CONFIG.declare('xisint', ConfigValue(
    #     default=None,
    #     domain=list_of_ints,
    #     description="List of inputs that should be treated as integers.",
    #     doc="Row vector of 0/1 flags that specify which input variables, if "
    #     "any, ALAMO should treat as integers. For integer inputs, ALAMO’s "
    #     "sampling will be restricted to integer values."))
    # CONFIG.declare('zisint', ConfigValue(
    #     default=None,
    #     domain=list_of_ints,
    #     description="List of outputs that should be treated as integers.",
    #     doc="Row vector of 0/1 flags that specify which output variables, if "
    #     "any, ALAMO should treat as integers. For integer variables, ALAMO’s "
    #     "model will include the rounding of a function to the nearest integer "
    #     "(equivalent to the nint function in Fortran.)"))

    # CONFIG.declare('tolrelmetric', ConfigValue(
    #     default=None,
    #     domain=list_of_floats,
    #     description="Relative tolerance for outputs.",
    #     doc="Relative convergence tolerance for the chosen fitness metric for "
    #     "the modeling of output variables. One per output variable, space "
    #     "separated. Incremental model building will stop if two consecutive "
    #     "iterations do not improve the chosen metric by at least this amount."
    #     ))
    # CONFIG.declare('tolabsmetric', ConfigValue(
    #     default=None,
    #     domain=list_of_floats,
    #     description="Absolute tolerance for outputs.",
    #     doc="Absolute convergence tolerance for the chosen fitness metric for "
    #     "the modeling of output variables. One per output variable, space "
    #     "separated. Incremental model building will stop if two consecutive "
    #     "iterations do not improve the chosen metric by at least this amount."
    #     ))
    # CONFIG.declare('tolmeanerror', ConfigValue(
    #     default=None,
    #     domain=list_of_floats,
    #     description="Convergence tolerance for mean errors in outputs.",
    #     doc="Row vector of convergence tolerances for mean errors in the "
    #     "modeling of output variables. One per output variable, space "
    #     "separated. Incremental model building will stop if tolmeanerror, "
    #     "tolrelmetric, or tolabsmetric is satisfied."))
    # CONFIG.declare('tolsse', ConfigValue(
    #     default=None,
    #     domain=float,
    #     description="Absolute tolerance on SSE",
    #     doc="Absolute tolerance on sum of square errors (SSE). ALAMO will "
    #     "terminate if it finds a solution whose SSE is within tolsse from "
    #     "the SSE of the full least squares problem."))

    # CONFIG.declare('mipoptca', ConfigValue(
    #     default=None,
    #     domain=float,
    #     description="Absolute tolerance for MIP."))
    # CONFIG.declare('mipoptcr', ConfigValue(
    #     default=None,
    #     domain=float,
    #     description="Relative tolerance for MIP."))
    # CONFIG.declare('linearerror', ConfigValue(
    #     default=None,
    #     domain=In([True, False]),
    #     description="If True, a linear objective is used when solving "
    #     "mixed-integer optimization problems; otherwise, a squared error will "
    #     "be employed."))
    # CONFIG.declare('GAMS', ConfigValue(
    #     default=None,
    #     domain=str,
    #     description="Complete path of GAMS executable (or name if GAMS is in "
    #     "the user path)."))
    # CONFIG.declare('GAMSSOLVER', ConfigValue(
    #     default=None,
    #     domain=str,
    #     description="Name of preferred GAMS solver for solving ALAMO’s "
    #     "mixed-integer quadratic subproblems. Special facilities have been "
    #     "implemented in ALAMO and BARON that make BARON the preferred "
    #     "selection for this option. However, any mixed-integer quadratic "
    #     "programming solver available under GAMS can be used."))
    # CONFIG.declare('solvemip', ConfigValue(
    #     default=None,
    #     domain=In([True, False]),
    #     description="Whether to use an optimizer to solve MIP."))
    # CONFIG.declare('log_output', ConfigValue(
    #     default=False,
    #     domain=In([True, False]),
    #     description="Whether to log ALAMO output."))