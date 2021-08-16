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
import os

from idaes.surrogate.alamopy_new import Alamopy, Modelers, Screener


@pytest.fixture
def alm_obj():
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
def test_get_path(alm_obj):
    alm_obj.config.temp_path = os.path.dirname(__file__)
    alm_obj.config.filename = "foo.alm"
    a, t = alm_obj.get_paths()
    assert a == "foo.alm"
    assert t == "foo.trc"
    assert os.getcwd() == os.path.dirname(__file__)


@pytest.mark.unit
def test_get_path_default(alm_obj):
    a, t = alm_obj.get_paths()
    assert a == "alamopy.alm"
    assert t == "alamopy.trc"
    assert os.getcwd() == os.path.dirname(__file__)[:-6]


@pytest.mark.unit
def test_writer_default(alm_obj):
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
def test_writer_trace(alm_obj):
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
def test_writer_vdata(alm_obj):
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
def test_writer_full_config(alm_obj):
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


@pytest.mark.unit
def test_read_trace(alm_obj):
    os.chdir("tests")
    trc = alm_obj.read_trace_file(trace_file="alamotrace.trc")

    mstring = (
        " z1 = 3.9999999999925446303450 * x1^2 - "
        "4.0000000000020765611453 * x2^2 - 2.0999999999859380039879 * x1^4 + "
        "4.0000000000043112180492 * x2^4 + 0.33333333332782633107172 * x1^6 + "
        "0.99999999999972988273811 * x1*x2")

    assert trc == {
        "#filename": "GUI_camel6.alm",
        "NINPUTS": "2",
        "NOUTPUTS": "1",
        "INITIALPOINTS": "10",
        "OUTPUT": "1",
        "SET": "0",
        "INITIALIZER": "3",
        "SAMPLER": "1",
        "MODELER": "1",
        "BUILDER": "1",
        "GREEDYBUILD": "T",
        "BACKSTEPPER": "0",
        "GREEDYBACK": "T",
        "REGULARIZER": "0",
        "SOLVEMIP": "F",
        "SSEOLR": "0.602E-28",
        "SSE": "0.976E-23",
        "RMSE": "0.988E-12",
        "R2": "1.00",
        "ModelSize": "6",
        "BIC": "-539.",
        "RIC": "32.5",
        "Cp": "2.00",
        "AICc": "-513.",
        "HQC": "-543.",
        "MSE": "0.325E-23",
        "SSEp": "0.976E-23",
        "MADp": "0.115E-07",
        "OLRTime": "0.46875000E-01",
        "numOLRs": "16384",
        "OLRoneCalls": "15",
        "OLRoneFails": "0",
        "OLRgsiCalls": "0",
        "OLRgsiFails": "0",
        "OLRdgelCalls": "16369",
        "OLRdgelFails": "0",
        "OLRclrCalls": "0",
        "OLRclrFails": "0",
        "OLRgmsCalls": "0",
        "OLRgmsFails": "0",
        "CLRTime": "0.0000000",
        "numCLRs": "0",
        "MIPTime": "0.0000000",
        "NumMIPs": "0",
        "LassoTime": "0.0000000",
        "Metric1Lasso": "0.17976931+309",
        "Metric2Lasso": "0.17976931+309",
        "LassoSuccess": "F",
        "LassoRed": "0.0000000",
        "nBasInitAct": "15",
        "nBas": "15",
        "SimTime": "0.0000000",
        "SimData": "0",
        "TotData": "10",
        "NdataConv": "0",
        "OtherTime": "0.46875000E-01",
        "NumIters": "1",
        "IterConv": "0",
        "TimeConv": "0.0000000",
        "Step0Time": "0.0000000",
        "Step1Time": "0.93750000E-01",
        "Step2Time": "0.0000000",
        "TotalTime": "0.93750000E-01",
        "AlamoStatus": "0",
        "AlamoVersion": "2021.5.8",
        "Model": mstring
        }
