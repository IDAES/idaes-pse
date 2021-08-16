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
        alm_obj.write_alm_to_stream(stream=stream)

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
        alm_obj.write_alm_to_stream(
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
        alm_obj.write_alm_to_stream(stream=stream)

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
        alm_obj.config.convpen = 42
        alm_obj.config.screener = Screener.SIS
        alm_obj.config.ncvf = 4
        alm_obj.config.sismult = 5

        alm_obj.config.maxiter = -50
        alm_obj.config.maxtime = 500
        alm_obj.config.datalimitterms = True
        alm_obj.config.maxterms = [-1, 2]
        alm_obj.config.minterms = [2, -1]
        alm_obj.config.numlimitbasis = False
        alm_obj.config.exclude = [1, 0]
        alm_obj.config.ignore = 1
        alm_obj.config.xisint = [1, 0]
        alm_obj.config.zisint = 0

        alm_obj.config.tolrelmetric = 8.1
        alm_obj.config.tolabsmetric = 9.2
        alm_obj.config.tolmeanerror = 10.3
        alm_obj.config.tolsse = 11.4

        alm_obj.config.mipoptca = 12.5
        alm_obj.config.mipoptcr = 13.6
        alm_obj.config.linearerror = False
        alm_obj.config.GAMS = "foo"
        alm_obj.config.GAMSSOLVER = "bar"
        alm_obj.config.solvemip = False

        stream = io.StringIO()
        alm_obj.write_alm_to_stream(stream=stream)

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
            "convpen 42.0\n"
            "screener 2\n"
            "ncvf 4\n"
            "sismult 5\n"
            "maxtime 500.0\n"
            "maxiter -50\n"
            "datalimitterms 1\n"
            "maxterms -1 2\n"
            "minterms 2 -1\n"
            "numlimitbasis 0\n"
            "exclude 1 0\n"
            "ignore 1\n"
            "xisint 1 0\n"
            "zisint 0\n"
            "tolrelmetric 8.1\n"
            "tolabsmetric 9.2\n"
            "tolmeanerror 10.3\n"
            "tolsse 11.4\n"
            "mipoptca 12.5\n"
            "mipoptcr 13.6\n"
            "linearerror 0\n"
            "GAMS foo\n"
            "GAMSSOLVER bar\n"
            "solvemip 0\n"
            "\nTRACE 1\n\n"
            "BEGIN_DATA\n"
            "1 5 10\n"
            "2 6 20\n"
            "3 7 30\n"
            "4 8 40\n"
            "END_DATA\n")
