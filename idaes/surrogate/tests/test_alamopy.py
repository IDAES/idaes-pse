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
from math import sin, cos, log, exp
from pathlib import Path
from io import StringIO

from pyomo.environ import Var, Constraint
from pyomo.common.tempfiles import TempfileManager

from idaes.surrogate.alamopy_new import \
    AlamoTrainer, AlamoObject, Modelers, Screener, alamo
from idaes.surrogate.surrogate_block import SurrogateBlock
from idaes.core.util.exceptions import ConfigurationError


dirpath = Path(__file__).parent.resolve()


# String representation of json output for testing
jstring = (
    '{"surrogate": {"z1": " z1 == 3.9999999999925446303450 * x1**2 - '
    '4.0000000000020765611453 * x2**2 - '
    '2.0999999999859380039879 * x1**4 + '
    '4.0000000000043112180492 * x2**4 + '
    '0.33333333332782633107172 * x1**6 + '
    '0.99999999999972988273811 * x1*x2"}, '
    '"input_labels": ["x1", "x2"], '
    '"output_labels": ["z1"], '
    '"input_bounds": {"x1": [0, 5], "x2": [0, 10]}}')


class TestAlamoTrainer:
    @pytest.fixture
    def alm_obj(self):
        alm_obj = AlamoTrainer(input_labels=["x1", "x2"], output_labels=["z1"])

        alm_obj.set_input_bounds({"x1": (0, 5), "x2": (0, 10)})

        alm_obj.set_training_data(np.array([[1, 5], [2, 6], [3, 7], [4, 8]]),
                                  np.array([[10], [20], [30], [40]]))

        return alm_obj

    @pytest.mark.unit
    def test_get_files(self, alm_obj):
        alm_obj.config.filename = "foo.alm"
        alm_obj.get_files()
        assert alm_obj._almfile == "foo.alm"
        assert alm_obj._trcfile == "foo.trc"

    @pytest.mark.unit
    def test_get_files_default(self, alm_obj):
        alm_obj.get_files()

        assert alm_obj._almfile is not None
        assert os.path.exists(alm_obj._almfile)
        assert str(alm_obj._trcfile).split(".")[0] == str(
            alm_obj._almfile).split(".")[0]

        alm_obj.remove_temp_files()

        # Check that we cleaned-up after ourselves
        assert not os.path.exists(alm_obj._almfile)

    @pytest.mark.unit
    def test_get_files_exists(self, alm_obj):
        alm_obj.config.filename = os.path.join(dirpath, "alamotrace.trc")

        with pytest.raises(
                FileExistsError,
                match="alamotrace.trc already exists. Either choose a new "
                "file name or set overwrite_files = True"):
            alm_obj.get_files()

    @pytest.mark.unit
    def test_writer_default(self, alm_obj):
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
            "NDATA 4\n\n"
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
    def test_writer_custom_basis(self, alm_obj):
        alm_obj.config.custom_basis_functions = [
            "sin(x1*x2)", "x1*tanh(x2)"]
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
            "NDATA 4\n\n"
            "NCUSTOMBAS 2\n"
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
            "END_DATA\n\n"
            "BEGIN_CUSTOMBAS\n"
            "sin(x1*x2)\n"
            "x1*tanh(x2)\n"
            "END_CUSTOMBAS\n")

    @pytest.mark.unit
    def test_writer_min_max_equal(self, alm_obj):
        alm_obj.set_input_bounds({"x1": (0, 5), "x2": (10, 10)})
        stream = io.StringIO()

        with pytest.raises(ConfigurationError,
                           match="ALAMO configuration error: upper and lower "
                           "bounds on input x2 are equal."):
            alm_obj.write_alm_to_stream(stream=stream)

    @pytest.mark.unit
    def test_writer_min_max_reversed(self, alm_obj):
        alm_obj.set_input_bounds({"x1": (0, 5), "x2": (15, 10)})
        stream = io.StringIO()

        with pytest.raises(ConfigurationError,
                           match="ALAMO configuration error: upper bound is "
                           "less than lower bound for input x2."):
            alm_obj.write_alm_to_stream(stream=stream)

    @pytest.mark.unit
    def test_writer_trace(self, alm_obj):
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
            "NDATA 4\n\n"
            "linfcns 1\n"
            "constant 1\n"
            "maxtime 1000.0\n"
            "numlimitbasis 1\n\n"
            "TRACE 1\n"
            "TRACEFNAME foo.bar\n\n"
            "BEGIN_DATA\n"
            "1 5 10\n"
            "2 6 20\n"
            "3 7 30\n"
            "4 8 40\n"
            "END_DATA\n")

    @pytest.mark.unit
    def test_writer_validation_data(self, alm_obj):
        alm_obj.set_validation_data(np.array([[2.5, 6.5]]),
                                    np.array([[25]], ndmin=2))

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
    def test_writer_full_config(self, alm_obj):
        alm_obj.config.xfactor = [1.2, 1.6]
        alm_obj.config.xscaling = True
        alm_obj.config.scalez = True

        alm_obj.config.monomialpower = [2]
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
        alm_obj.config.print_to_screen = True

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
            "NDATA 4\n\n"
            "MONO 1\n"
            "MULTI2 3\n"
            "MULTI3 4\n"
            "RATIOS 5\n"
            "xfactor 1.2 1.6\n"
            "xscaling 1\n"
            "scalez 1\n"
            "monomialpower 2\n"
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
            "print_to_screen 1\n"
            "\nTRACE 1\n\n"
            "BEGIN_DATA\n"
            "1 5 10\n"
            "2 6 20\n"
            "3 7 30\n"
            "4 8 40\n"
            "END_DATA\n")

    @pytest.mark.component
    def test_file_writer(self, alm_obj):
        alm_obj.get_files()
        alm_obj.write_alm_file()

        with open(alm_obj._almfile, "r") as f:
            fcont = f.read()
        f.close()
        alm_obj.remove_temp_files()

        assert fcont == (
            f"# IDAES Alamopy input file\n"
            f"NINPUTS 2\n"
            f"NOUTPUTS 1\n"
            f"XLABELS x1 x2\n"
            f"ZLABELS z1\n"
            f"XMIN 0 0\n"
            f"XMAX 5 10\n"
            f"NDATA 4\n\n"
            f"linfcns 1\n"
            f"constant 1\n"
            f"maxtime 1000.0\n"
            f"numlimitbasis 1\n\n"
            f"TRACE 1\n"
            f"TRACEFNAME {alm_obj._trcfile}\n\n"
            f"BEGIN_DATA\n"
            f"1 5 10\n"
            f"2 6 20\n"
            f"3 7 30\n"
            f"4 8 40\n"
            f"END_DATA\n")

        # Check for clean up
        assert not os.path.exists(alm_obj._almfile)

    @pytest.mark.unit
    @pytest.mark.skipif(not alamo.available(), reason="ALAMO not available")
    def test_call_alamo_empty(self, alm_obj):
        alm_obj._almfile = "test"
        with pytest.raises(RuntimeError,
                           match="ALAMO executable returned non-zero return "
                           "code. Check the ALAMO output for more "
                           "information."):
            alm_obj.call_alamo()

        assert alm_obj._temp_context is None

    @pytest.mark.component
    @pytest.mark.skipif(not alamo.available(), reason="ALAMO not available")
    def test_call_alamo_w_input(self, alm_obj):
        cwd = os.getcwd()

        alm_obj._temp_context = TempfileManager.new_context()
        alm_obj._almfile = os.path.join(dirpath, "alamo_test.alm")
        alm_obj._wrkdir = dirpath

        alm_obj._temp_context.add_tempfile(
            os.path.join(dirpath, "GUI_trace.trc"), exists=False)
        alm_obj._temp_context.add_tempfile(
            os.path.join(dirpath, "alamo_test.lst"), exists=False)
        rc, almlog = alm_obj.call_alamo()
        assert rc == 0
        assert "Normal termination" in almlog

        # Check for clean up
        alm_obj._temp_context.release(remove=True)
        assert not os.path.exists(os.path.join(dirpath, "GUI_trace.trc"))
        assert not os.path.exists(os.path.join(dirpath, "alamo_test.lst"))
        assert cwd == os.getcwd()

    @pytest.mark.unit
    def test_read_trace_single(self, alm_obj):
        alm_obj._trcfile = os.path.join(dirpath, "alamotrace.trc")
        trc = alm_obj.read_trace_file()

        mdict = {'z1': (' z1 == 3.9999999999925446303450 * x1**2 - '
                        '4.0000000000020765611453 * x2**2 - '
                        '2.0999999999859380039879 * x1**4 + '
                        '4.0000000000043112180492 * x2**4 + '
                        '0.33333333332782633107172 * x1**6 + '
                        '0.99999999999972988273811 * x1*x2')}

        assert trc == {
            'filename': 'GUI_camel6.alm',
            'NINPUTS': '2',
            'NOUTPUTS': '1',
            'INITIALPOINTS': '10',
            'OUTPUT': {'z1': '1'},
            'SET': '0',
            'INITIALIZER': '3',
            'SAMPLER': '1',
            'MODELER': '1',
            'BUILDER': '1',
            'GREEDYBUILD': 'T',
            'BACKSTEPPER': '0',
            'GREEDYBACK': 'T',
            'REGULARIZER': '0',
            'SOLVEMIP': 'F',
            'SSEOLR': {'z1': '0.602E-28'},
            'SSE': {'z1': '0.976E-23'},
            'RMSE': {'z1': '0.988E-12'},
            'R2': {'z1': '1.00'},
            'ModelSize': {'z1': '6'},
            'BIC': {'z1': '-539.'},
            'RIC': {'z1': '32.5'},
            'Cp': {'z1': '2.00'},
            'AICc': {'z1': '-513.'},
            'HQC': {'z1': '-543.'},
            'MSE': {'z1': '0.325E-23'},
            'SSEp': {'z1': '0.976E-23'},
            'MADp': {'z1': '0.115E-07'},
            'OLRTime': {'z1': '0.46875000E-01'},
            'numOLRs': {'z1': '16384'},
            'OLRoneCalls': {'z1': '15'},
            'OLRoneFails': {'z1': '0'},
            'OLRgsiCalls': {'z1': '0'},
            'OLRgsiFails': {'z1': '0'},
            'OLRdgelCalls': {'z1': '16369'},
            'OLRdgelFails': {'z1': '0'},
            'OLRclrCalls': {'z1': '0'},
            'OLRclrFails': {'z1': '0'},
            'OLRgmsCalls': {'z1': '0'},
            'OLRgmsFails': {'z1': '0'},
            'CLRTime': {'z1': '0.0000000'},
            'numCLRs': {'z1': '0'},
            'MIPTime': {'z1': '0.0000000'},
            'NumMIPs': {'z1': '0'},
            'LassoTime': {'z1': '0.0000000'},
            'Metric1Lasso': {'z1': '0.17976931+309'},
            'Metric2Lasso': {'z1': '0.17976931+309'},
            'LassoSuccess': {'z1': 'F'},
            'LassoRed': {'z1': '0.0000000'},
            'nBasInitAct': {'z1': '15'},
            'nBas': {'z1': '15'},
            'SimTime': {'z1': '0.0000000'},
            'SimData': {'z1': '0'},
            'TotData': {'z1': '10'},
            'NdataConv': {'z1': '0'},
            'OtherTime': {'z1': '0.46875000E-01'},
            'NumIters': {'z1': '1'},
            'IterConv': {'z1': '0'},
            'TimeConv': {'z1': '0.0000000'},
            'Step0Time': {'z1': '0.0000000'},
            'Step1Time': {'z1': '0.93750000E-01'},
            'Step2Time': {'z1': '0.0000000'},
            'TotalTime': {'z1': '0.93750000E-01'},
            'AlamoStatus': {'z1': '0'},
            'AlamoVersion': {'z1': '2021.5.8'},
            'Model': mdict}

    @pytest.mark.unit
    def test_read_trace_multi(self, alm_obj):
        alm_obj._output_labels = ["z1", "z2"]

        alm_obj._trcfile = os.path.join(dirpath, "alamotrace2.trc")
        trc = alm_obj.read_trace_file()

        mdict = {
            'z1': (' z1 == 3.9999999999925446303450 * x1**2 - '
                   '4.0000000000020765611453 * x2**2 - '
                   '2.0999999999859380039879 * x1**4 + '
                   '4.0000000000043112180492 * x2**4 + '
                   '0.33333333332782633107172 * x1**6 + '
                   '0.99999999999972988273811 * x1*x2'),
            'z2': (' z2 == 0.72267799849202937756409E-001 * x1 + '
                   '0.68451684753912708791823E-001 * x2 + '
                   '1.0677896915911471165117 * x1**2 - '
                   '0.70576464806224348258468 * x2**2 - '
                   '0.40286283566554989543640E-001 * x1**3 + '
                   '0.67785668021684807038607E-002 * x2**5 - '
                   '0.14017881460354553180281 * x1**6 + '
                   '0.77207049441576647286212 * x2**6 + '
                   '0.42143309951518070910481 * x1*x2 - '
                   '0.41818729807213093907503E-001')}

        assert trc == {
            'filename': 'GUI_camel6_twooutputs.alm',
            'NINPUTS': '2',
            'NOUTPUTS': '2',
            'INITIALPOINTS': '10',
            'OUTPUT': {'z1': '1', 'z2': '2'},
            'SET': '0',
            'INITIALIZER': '3',
            'SAMPLER': '1',
            'MODELER': '1',
            'BUILDER': '1',
            'GREEDYBUILD': 'T',
            'BACKSTEPPER': '0',
            'GREEDYBACK': 'T',
            'REGULARIZER': '0',
            'SOLVEMIP': 'F',
            'SSEOLR': {'z1': '0.602E-28', 'z2': '0.280E-30'},
            'SSE': {'z1': '0.976E-23', 'z2': '0.280E-30'},
            'RMSE': {'z1': '0.988E-12', 'z2': '0.167E-15'},
            'R2': {'z1': '1.00', 'z2': '1.00'},
            'ModelSize': {'z1': '6', 'z2': '10'},
            'BIC': {'z1': '-539.', 'z2': '-704.'},
            'RIC': {'z1': '32.5', 'z2': '54.2'},
            'Cp': {'z1': '2.00', 'z2': '10.0'},
            'AICc': {'z1': '-513.', 'z2': '-707.'},
            'HQC': {'z1': '-543.', 'z2': '-710.'},
            'MSE': {'z1': '0.325E-23', 'z2': '0.280E-30'},
            'SSEp': {'z1': '0.976E-23', 'z2': '0.280E-30'},
            'MADp': {'z1': '0.115E-07', 'z2': '0.118E-11'},
            'OLRTime': {'z1': '0.15625000E-01', 'z2': '0.93750000E-01'},
            'numOLRs': {'z1': '16384', 'z2': '30827'},
            'OLRoneCalls': {'z1': '30', 'z2': '30'},
            'OLRoneFails': {'z1': '0', 'z2': '0'},
            'OLRgsiCalls': {'z1': '0', 'z2': '0'},
            'OLRgsiFails': {'z1': '0', 'z2': '0'},
            'OLRdgelCalls': {'z1': '47181', 'z2': '47181'},
            'OLRdgelFails': {'z1': '0', 'z2': '0'},
            'OLRclrCalls': {'z1': '0', 'z2': '0'},
            'OLRclrFails': {'z1': '0', 'z2': '0'},
            'OLRgmsCalls': {'z1': '0', 'z2': '0'},
            'OLRgmsFails': {'z1': '0', 'z2': '0'},
            'CLRTime': {'z1': '0.0000000', 'z2': '0.0000000'},
            'numCLRs': {'z1': '0', 'z2': '0'},
            'MIPTime': {'z1': '0.0000000', 'z2': '0.0000000'},
            'NumMIPs': {'z1': '0', 'z2': '0'},
            'LassoTime': {'z1': '0.0000000', 'z2': '0.0000000'},
            'Metric1Lasso': {'z1': '0.17976931+309', 'z2': '0.17976931+309'},
            'Metric2Lasso': {'z1': '0.17976931+309', 'z2': '0.17976931+309'},
            'LassoSuccess': {'z1': 'F', 'z2': 'F'},
            'LassoRed': {'z1': '0.0000000', 'z2': '0.0000000'},
            'nBasInitAct': {'z1': '15', 'z2': '15'},
            'nBas': {'z1': '15', 'z2': '15'},
            'SimTime': {'z1': '0.0000000', 'z2': '0.0000000'},
            'SimData': {'z1': '0', 'z2': '0'},
            'TotData': {'z1': '10', 'z2': '10'},
            'NdataConv': {'z1': '0', 'z2': '0'},
            'OtherTime': {'z1': '0.31250000E-01', 'z2': '0.31250000E-01'},
            'NumIters': {'z1': '1', 'z2': '1'},
            'IterConv': {'z1': '0', 'z2': '0'},
            'TimeConv': {'z1': '0.0000000', 'z2': '0.0000000'},
            'Step0Time': {'z1': '0.0000000', 'z2': '0.0000000'},
            'Step1Time': {'z1': '0.15625000', 'z2': '0.15625000'},
            'Step2Time': {'z1': '0.0000000', 'z2': '0.0000000'},
            'TotalTime': {'z1': '0.46875000E-01', 'z2': '0.10937500'},
            'AlamoStatus': {'z1': '0', 'z2': '0'},
            'AlamoVersion': {'z1': '2021.5.8', 'z2': '2021.5.8'},
            'Model': mdict}

    @pytest.mark.unit
    def test_read_trace_number_mismatch(self, alm_obj):
        alm_obj._trcfile = os.path.join(dirpath, "alamotrace2.trc")
        with pytest.raises(RuntimeError,
                           match="Mismatch when reading ALAMO trace file. "
                           "Expected OUTPUT = 1, found 2."):
            alm_obj.read_trace_file()

    @pytest.mark.unit
    def test_read_trace_label_mismatch(self, alm_obj):
        alm_obj._output_labels = ["z1", "z3"]

        alm_obj._trcfile = os.path.join(dirpath, "alamotrace2.trc")
        with pytest.raises(RuntimeError,
                           match="Mismatch when reading ALAMO trace file. "
                           "Label of output variable in expression "
                           "\(z2\) does not match expected label \(z3\)."):
            alm_obj.read_trace_file()

    @pytest.mark.unit
    def test_populate_results(self, alm_obj):
        alm_obj._trcfile = os.path.join(dirpath, "alamotrace.trc")
        trc = alm_obj.read_trace_file()
        alm_obj.populate_results(trc)

        mdict = {'z1': (' z1 == 3.9999999999925446303450 * x1**2 - '
                        '4.0000000000020765611453 * x2**2 - '
                        '2.0999999999859380039879 * x1**4 + '
                        '4.0000000000043112180492 * x2**4 + '
                        '0.33333333332782633107172 * x1**6 + '
                        '0.99999999999972988273811 * x1*x2')}

        assert alm_obj._results == {
            'filename': 'GUI_camel6.alm',
            'NINPUTS': '2',
            'NOUTPUTS': '1',
            'INITIALPOINTS': '10',
            'OUTPUT': {'z1': '1'},
            'SET': '0',
            'INITIALIZER': '3',
            'SAMPLER': '1',
            'MODELER': '1',
            'BUILDER': '1',
            'GREEDYBUILD': 'T',
            'BACKSTEPPER': '0',
            'GREEDYBACK': 'T',
            'REGULARIZER': '0',
            'SOLVEMIP': 'F',
            'SSEOLR': {'z1': '0.602E-28'},
            'SSE': {'z1': '0.976E-23'},
            'RMSE': {'z1': '0.988E-12'},
            'R2': {'z1': '1.00'},
            'ModelSize': {'z1': '6'},
            'BIC': {'z1': '-539.'},
            'RIC': {'z1': '32.5'},
            'Cp': {'z1': '2.00'},
            'AICc': {'z1': '-513.'},
            'HQC': {'z1': '-543.'},
            'MSE': {'z1': '0.325E-23'},
            'SSEp': {'z1': '0.976E-23'},
            'MADp': {'z1': '0.115E-07'},
            'OLRTime': {'z1': '0.46875000E-01'},
            'numOLRs': {'z1': '16384'},
            'OLRoneCalls': {'z1': '15'},
            'OLRoneFails': {'z1': '0'},
            'OLRgsiCalls': {'z1': '0'},
            'OLRgsiFails': {'z1': '0'},
            'OLRdgelCalls': {'z1': '16369'},
            'OLRdgelFails': {'z1': '0'},
            'OLRclrCalls': {'z1': '0'},
            'OLRclrFails': {'z1': '0'},
            'OLRgmsCalls': {'z1': '0'},
            'OLRgmsFails': {'z1': '0'},
            'CLRTime': {'z1': '0.0000000'},
            'numCLRs': {'z1': '0'},
            'MIPTime': {'z1': '0.0000000'},
            'NumMIPs': {'z1': '0'},
            'LassoTime': {'z1': '0.0000000'},
            'Metric1Lasso': {'z1': '0.17976931+309'},
            'Metric2Lasso': {'z1': '0.17976931+309'},
            'LassoSuccess': {'z1': 'F'},
            'LassoRed': {'z1': '0.0000000'},
            'nBasInitAct': {'z1': '15'},
            'nBas': {'z1': '15'},
            'SimTime': {'z1': '0.0000000'},
            'SimData': {'z1': '0'},
            'TotData': {'z1': '10'},
            'NdataConv': {'z1': '0'},
            'OtherTime': {'z1': '0.46875000E-01'},
            'NumIters': {'z1': '1'},
            'IterConv': {'z1': '0'},
            'TimeConv': {'z1': '0.0000000'},
            'Step0Time': {'z1': '0.0000000'},
            'Step1Time': {'z1': '0.93750000E-01'},
            'Step2Time': {'z1': '0.0000000'},
            'TotalTime': {'z1': '0.93750000E-01'},
            'AlamoStatus': {'z1': '0'},
            'AlamoVersion': {'z1': '2021.5.8'},
            'Model': mdict}

    @pytest.mark.unit
    def test_build_surrogate_object(self, alm_obj):
        alm_obj._trcfile = os.path.join(dirpath, "alamotrace.trc")
        trc = alm_obj.read_trace_file()
        alm_obj.populate_results(trc)
        alm_obj.build_surrogate_object()

        assert isinstance(alm_obj._surrogate, AlamoObject)
        assert alm_obj._surrogate._surrogate == {
            'z1': (' z1 == 3.9999999999925446303450 * x1**2 - '
                   '4.0000000000020765611453 * x2**2 - '
                   '2.0999999999859380039879 * x1**4 + '
                   '4.0000000000043112180492 * x2**4 + '
                   '0.33333333332782633107172 * x1**6 + '
                   '0.99999999999972988273811 * x1*x2')}
        assert alm_obj._surrogate._input_labels == ["x1", "x2"]
        assert alm_obj._surrogate._output_labels == ["z1"]
        assert alm_obj._surrogate._input_bounds == {
            "x1": (0, 5), "x2": (0, 10)}


class TestAlamoObject():
    @pytest.fixture
    def alm_surr1(self):
        surrogate = {
            'z1': (' z1 == 3.9999999999925446303450 * x1**2 - '
                   '4.0000000000020765611453 * x2**2 - '
                   '2.0999999999859380039879 * x1**4 + '
                   '4.0000000000043112180492 * x2**4 + '
                   '0.33333333332782633107172 * x1**6 + '
                   '0.99999999999972988273811 * x1*x2')}
        input_labels = ["x1", "x2"]
        output_labels = ["z1"]
        input_bounds = {"x1": (0, 5), "x2": (0, 10)}

        alm_surr1 = AlamoObject(
            surrogate, input_labels, output_labels, input_bounds)

        return alm_surr1

    @pytest.mark.unit
    def test_evaluate_surrogate(self, alm_surr1):
        x = [-2, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0,
             0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]

        inputs = np.array([np.tile(x, len(x)), np.repeat(x, len(x))])

        out = alm_surr1.evaluate_surrogate(inputs)
        for i in range(inputs.shape[1]):
            assert pytest.approx(out[0, i], rel=1e-8) == (
                3.9999999999925446303450 * inputs[0, i]**2 -
                4.0000000000020765611453 * inputs[1, i]**2 -
                2.0999999999859380039879 * inputs[0, i]**4 +
                4.0000000000043112180492 * inputs[1, i]**4 +
                0.33333333332782633107172 * inputs[0, i]**6 +
                0.99999999999972988273811 * inputs[0, i]*inputs[1, i])

    @pytest.mark.unit
    def test_populate_block(self, alm_surr1):
        blk = SurrogateBlock(concrete=True)

        blk.build_model(alm_surr1)
        blk.display()

        assert isinstance(blk._inputs, Var)
        assert blk._inputs["x1"].bounds == (0, 5)
        assert blk._inputs["x2"].bounds == (0, 10)
        assert isinstance(blk._outputs, Var)
        assert blk._outputs["z1"].bounds == (None, None)
        assert isinstance(blk.alamo_constraint, Constraint)
        assert len(blk.alamo_constraint) == 1
        assert str(blk.alamo_constraint["z1"].body) == (
            "_outputs[z1] - ("
            "3.9999999999925446*_inputs[x1]**2 - "
            "4.000000000002077*_inputs[x2]**2 - "
            "2.099999999985938*_inputs[x1]**4 + "
            "4.000000000004311*_inputs[x2]**4 + "
            "0.33333333332782633*_inputs[x1]**6 + "
            "0.9999999999997299*_inputs[x1]*_inputs[x2])")

    # @pytest.mark.unit
    # def test_populate_block_indexed(self, alm_surr1):
    #     blk = Block(concrete=True)

    #     alm_surr1.populate_block(blk, index_set=[1, 2, 3])

    #     assert isinstance(blk.x1, Var)
    #     assert blk.x1.is_indexed()
    #     assert len(blk.x1) == 3
    #     assert isinstance(blk.x2, Var)
    #     assert blk.x2.is_indexed()
    #     assert len(blk.x2) == 3
    #     assert isinstance(blk.z1, Var)
    #     assert blk.z1.is_indexed()
    #     assert len(blk.z1) == 3
    #     assert isinstance(blk.alamo_constraint, Constraint)
    #     assert blk.alamo_constraint.is_indexed()
    #     assert len(blk.alamo_constraint) == 3

    #     for i in [1, 2, 3]:
    #         assert blk.x1[i].bounds == (0, 5)
    #         assert blk.x2[i].bounds == (0, 10)
    #         assert blk.z1[i].bounds == (None, None)
    #         assert str(blk.alamo_constraint["z1", i].body) == (
    #             f"z1[{i}] - (3.9999999999925446*x1[{i}]**2 - "
    #             f"4.000000000002077*x2[{i}]**2 - "
    #             f"2.099999999985938*x1[{i}]**4 + "
    #             f"4.000000000004311*x2[{i}]**4 + "
    #             f"0.33333333332782633*x1[{i}]**6 + "
    #             f"0.9999999999997299*x1[{i}]*x2[{i}])")

    @pytest.fixture
    def alm_surr2(self):
        surrogate = {
            'z1': ('z1 == 3.9999999999925446303450 * x1**2 - '
                   '4.0000000000020765611453 * x2**2 - '
                   '2.0999999999859380039879 * x1**4 + '
                   '4.0000000000043112180492 * x2**4 + '
                   '0.33333333332782633107172 * x1**6 + '
                   '0.99999999999972988273811 * x1*x2'),
            'z2': ('z2 == 0.72267799849202937756409E-001 * x1 + '
                   '0.68451684753912708791823E-001 * x2 + '
                   '1.0677896915911471165117 * x1**2 - '
                   '0.70576464806224348258468 * x2**2 - '
                   '0.40286283566554989543640E-001 * x1**3 + '
                   '0.67785668021684807038607E-002 * x2**5 - '
                   '0.14017881460354553180281 * x1**6 + '
                   '0.77207049441576647286212 * x2**6 + '
                   '0.42143309951518070910481 * x1*x2 - '
                   '0.41818729807213093907503E-001')}
        input_labels = ["x1", "x2"]
        output_labels = ["z1", "z2"]

        alm_surr2 = AlamoObject(surrogate, input_labels, output_labels)

        return alm_surr2

    @pytest.mark.unit
    def test_evaluate_surrogate_multi(self, alm_surr2):
        x = [-2, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0,
             0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]

        inputs = np.array([np.tile(x, len(x)), np.repeat(x, len(x))])

        out = alm_surr2.evaluate_surrogate(inputs)
        for i in range(inputs.shape[1]):
            assert pytest.approx(out[0, i], rel=1e-8) == (
                3.9999999999925446303450 * inputs[0, i]**2 -
                4.0000000000020765611453 * inputs[1, i]**2 -
                2.0999999999859380039879 * inputs[0, i]**4 +
                4.0000000000043112180492 * inputs[1, i]**4 +
                0.33333333332782633107172 * inputs[0, i]**6 +
                0.99999999999972988273811 * inputs[0, i]*inputs[1, i])
            assert pytest.approx(out[1, i], rel=1e-8) == (
                0.72267799849202937756409E-001 * inputs[0, i] +
                0.68451684753912708791823E-001 * inputs[1, i] +
                1.0677896915911471165117 * inputs[0, i]**2 -
                0.70576464806224348258468 * inputs[1, i]**2 -
                0.40286283566554989543640E-001 * inputs[0, i]**3 +
                0.67785668021684807038607E-002 * inputs[1, i]**5 -
                0.14017881460354553180281 * inputs[0, i]**6 +
                0.77207049441576647286212 * inputs[1, i]**6 +
                0.42143309951518070910481 * inputs[0, i]*inputs[1, i] -
                0.41818729807213093907503E-001)

    @pytest.mark.unit
    def test_populate_block_multi(self, alm_surr2):
        blk = SurrogateBlock(concrete=True)

        blk.build_model(alm_surr2)

        assert isinstance(blk._inputs, Var)
        assert blk._inputs["x1"].bounds == (None, None)
        assert blk._inputs["x2"].bounds == (None, None)
        assert isinstance(blk._outputs, Var)
        assert blk._outputs["z1"].bounds == (None, None)
        assert blk._outputs["z2"].bounds == (None, None)
        assert isinstance(blk.alamo_constraint, Constraint)
        assert len(blk.alamo_constraint) == 2
        assert str(blk.alamo_constraint["z1"].body) == (
            "_outputs[z1] - ("
            "3.9999999999925446*_inputs[x1]**2 - "
            "4.000000000002077*_inputs[x2]**2 - "
            "2.099999999985938*_inputs[x1]**4 + "
            "4.000000000004311*_inputs[x2]**4 + "
            "0.33333333332782633*_inputs[x1]**6 + "
            "0.9999999999997299*_inputs[x1]*_inputs[x2])")
        assert str(blk.alamo_constraint["z2"].body) == (
            "_outputs[z2] - ("
            "0.07226779984920294*_inputs[x1] + "
            "0.06845168475391271*_inputs[x2] + "
            "1.0677896915911471*_inputs[x1]**2 - "
            "0.7057646480622435*_inputs[x2]**2 - "
            "0.04028628356655499*_inputs[x1]**3 + "
            "0.006778566802168481*_inputs[x2]**5 - "
            "0.14017881460354553*_inputs[x1]**6 + "
            "0.7720704944157665*_inputs[x2]**6 + "
            "0.4214330995151807*_inputs[x1]*_inputs[x2] - "
            "0.041818729807213094)")

    # @pytest.mark.unit
    # def test_populate_block_multi_indexed(self, alm_surr2):
    #     blk = Block(concrete=True)

    #     alm_surr2.populate_block(blk, index_set=[1, 2, 3])

    #     assert isinstance(blk.x1, Var)
    #     assert blk.x1.is_indexed()
    #     assert len(blk.x1) == 3
    #     assert isinstance(blk.x2, Var)
    #     assert blk.x2.is_indexed()
    #     assert len(blk.x2) == 3
    #     assert isinstance(blk.z1, Var)
    #     assert blk.z1.is_indexed()
    #     assert len(blk.z1) == 3
    #     assert isinstance(blk.alamo_constraint, Constraint)
    #     assert blk.alamo_constraint.is_indexed()
    #     assert len(blk.alamo_constraint) == 6

    #     for i in [1, 2, 3]:
    #         assert blk.x1[i].bounds == (None, None)
    #         assert blk.x2[i].bounds == (None, None)
    #         assert blk.z1[i].bounds == (None, None)
    #         assert str(blk.alamo_constraint["z1", i].body) == (
    #             f"z1[{i}] - (3.9999999999925446*x1[{i}]**2 - "
    #             f"4.000000000002077*x2[{i}]**2 - "
    #             f"2.099999999985938*x1[{i}]**4 + "
    #             f"4.000000000004311*x2[{i}]**4 + "
    #             f"0.33333333332782633*x1[{i}]**6 + "
    #             f"0.9999999999997299*x1[{i}]*x2[{i}])")
    #         assert str(blk.alamo_constraint["z2", i].body) == (
    #             f"z2[{i}] - (0.07226779984920294*x1[{i}] + "
    #             f"0.06845168475391271*x2[{i}] + "
    #             f"1.0677896915911471*x1[{i}]**2 - "
    #             f"0.7057646480622435*x2[{i}]**2 - "
    #             f"0.04028628356655499*x1[{i}]**3 + "
    #             f"0.006778566802168481*x2[{i}]**5 - "
    #             f"0.14017881460354553*x1[{i}]**6 + "
    #             f"0.7720704944157665*x2[{i}]**6 + "
    #             f"0.4214330995151807*x1[{i}]*x2[{i}] - "
    #             f"0.041818729807213094)")

    @pytest.fixture
    def alm_surr3(self):
        surrogate = {
            'z1': (' z1 == 2*sin(x1**2) - 3*cos(x2**3) - '
                   '4*log(x1**4) + 5*exp(x2**5)')}
        input_labels = ["x1", "x2"]
        output_labels = ["z1"]

        alm_surr3 = AlamoObject(surrogate, input_labels, output_labels)

        return alm_surr3

    @pytest.mark.unit
    def test_evaluate_surrogate_funcs(self, alm_surr3):
        x = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]

        inputs = np.array([np.tile(x, len(x)), np.repeat(x, len(x))])

        out = alm_surr3.evaluate_surrogate(inputs)
        for i in range(inputs.shape[1]):
            assert pytest.approx(out[0, i], rel=1e-8) == (
                2*sin(inputs[0, i]**2) - 3*cos(inputs[1, i]**3) -
                4*log(inputs[0, i]**4) + 5*exp(inputs[1, i]**5))

    @pytest.mark.unit
    def test_populate_block_funcs(self, alm_surr3):
        blk = SurrogateBlock(concrete=True)

        blk.build_model(alm_surr3)

        assert isinstance(blk._inputs, Var)
        assert isinstance(blk._outputs, Var)
        assert isinstance(blk.alamo_constraint, Constraint)
        assert len(blk.alamo_constraint) == 1
        assert str(blk.alamo_constraint["z1"].body) == (
            "_outputs[z1] - (2*sin(_inputs[x1]**2) - 3*cos(_inputs[x2]**3) - "
            "4*log(_inputs[x1]**4) + 5*exp(_inputs[x2]**5))")

    @pytest.mark.unit
    def test_to_json(self, alm_surr1):
        stream = StringIO()
        alm_surr1.to_json(stream)
        assert stream.getvalue() == jstring

    @pytest.mark.unit
    def test_from_json(self):
        alm_surr = AlamoObject({}, [], [])

        alm_surr.from_json(jstring)

        assert alm_surr._surrogate == {
            "z1": ' z1 == 3.9999999999925446303450 * x1**2 - '
            '4.0000000000020765611453 * x2**2 - '
            '2.0999999999859380039879 * x1**4 + '
            '4.0000000000043112180492 * x2**4 + '
            '0.33333333332782633107172 * x1**6 + '
            '0.99999999999972988273811 * x1*x2'}
        assert alm_surr._input_labels == ["x1", "x2"]
        assert alm_surr._output_labels == ["z1"]
        assert alm_surr._input_bounds == {"x1": (0, 5), "x2": (0, 10)}

    @pytest.mark.unit
    def test_save_load(self, alm_surr1):
        with TempfileManager as tf:
            fname = tf.create_tempfile(suffix=".json")
            alm_surr1.save(fname, overwrite=True)

            assert os.path.isfile(fname)

            with open(fname, "r") as f:
                js = f.read()
            f.close()

            alm_load = AlamoObject.load(fname)

        # Check file contents
        assert js == jstring

        # Check loaded object
        assert isinstance(alm_load, AlamoObject)
        assert alm_load._surrogate == {
            "z1": ' z1 == 3.9999999999925446303450 * x1**2 - '
            '4.0000000000020765611453 * x2**2 - '
            '2.0999999999859380039879 * x1**4 + '
            '4.0000000000043112180492 * x2**4 + '
            '0.33333333332782633107172 * x1**6 + '
            '0.99999999999972988273811 * x1*x2'}
        assert alm_load._input_labels == ["x1", "x2"]
        assert alm_load._output_labels == ["z1"]
        assert alm_load._input_bounds == {"x1": (0, 5), "x2": (0, 10)}

        # Check for clean up
        assert not os.path.isfile(fname)

    @pytest.mark.unit
    def test_save_no_overwrite(self, alm_surr1):
        with TempfileManager as tf:
            fname = tf.create_tempfile(suffix=".json")

            with pytest.raises(FileExistsError):
                alm_surr1.save(fname)

        # Check for clean up
        assert not os.path.isfile(fname)


@pytest.mark.integration
def test_workflow():
    # Test end-to-end workflow with a simple problem.
    alm_obj = AlamoTrainer(input_labels=["x1", "x2"], output_labels=["z1"])

    alm_obj.set_input_bounds({"x1": (-1.5, 1.5), "x2": (-1.5, 1.5)})

    alm_obj.set_training_data(
        np.array([[0.353837234435, 0.99275270941666],
                  [0.904978848612, -0.746908518721],
                  [0.643706630938, -0.617496599522],
                  [1.29881420688, 0.305594881575],
                  [1.35791650867, 0.351045058258],
                  [0.938369314089, -0.525167416293],
                  [-1.46593541641, 0.383902178482],
                  [-0.374378293218, -0.689730440659],
                  [0.690326213554, 0.569364994374],
                  [-0.961163301329, 0.499471920546]]),
        np.array([[0.762878272854], [0.387963718723], [-0.0205375902284],
                  [2.43011137696], [2.36989368612], [0.829756159423],
                  [1.14054797964], [-0.219122783909], [0.982068847698],
                  [0.936855365038]],
                 ndmin=2))

    alm_obj.config.linfcns = True
    alm_obj.config.monomialpower = [2, 3, 4, 5, 6]
    alm_obj.config.multi2power = [1, 2]

    rc, almlog = alm_obj.train_surrogate()

    # Check execution
    assert rc == 0
    assert "Normal termination" in almlog

    # Check temp file clean up
    assert alm_obj._temp_context is None
    assert not os.path.exists(alm_obj._almfile)
    assert not os.path.exists(alm_obj._trcfile)

    # Check results
    assert alm_obj._results is not None
    assert alm_obj._results['NINPUTS'] == '2'
    assert alm_obj._results['NOUTPUTS'] == '1'
    assert alm_obj._results['SSEOLR'] == {'z1': '0.373E-29'}
    assert alm_obj._results['SSE'] == {'z1': '0.976E-23'}
    assert alm_obj._results['RMSE'] == {'z1': '0.988E-12'}
    assert alm_obj._results['R2'] == {'z1': '1.00'}
    assert alm_obj._results['ModelSize'] == {'z1': '6'}
    assert alm_obj._results['BIC'] == {'z1': '-539.'}
    assert alm_obj._results['RIC'] == {'z1': '32.5'}
    assert alm_obj._results['Cp'] == {'z1': '2.00'}
    assert alm_obj._results['AICc'] == {'z1': '-513.'}
    assert alm_obj._results['HQC'] == {'z1': '-543.'}
    assert alm_obj._results['MSE'] == {'z1': '0.325E-23'}
    assert alm_obj._results['SSEp'] == {'z1': '0.976E-23'}
    assert alm_obj._results['MADp'] == {'z1': '0.115E-07'}

    assert alm_obj._results['Model'] == {
            "z1": " z1 == 3.9999999999925432980774 * x1**2 - "
            "4.0000000000020792256805 * x2**2 - "
            "2.0999999999859380039879 * x1**4 + "
            "4.0000000000043085535140 * x2**4 + "
            "0.33333333332782683067208 * x1**6 + "
            "0.99999999999973088193883 * x1*x2"}

    assert isinstance(alm_obj._surrogate, AlamoObject)
    assert alm_obj._surrogate._surrogate == {
        'z1': ' z1 == 3.9999999999925432980774 * x1**2 - '
        '4.0000000000020792256805 * x2**2 - '
        '2.0999999999859380039879 * x1**4 + '
        '4.0000000000043085535140 * x2**4 + '
        '0.33333333332782683067208 * x1**6 + '
        '0.99999999999973088193883 * x1*x2'}
    assert alm_obj._surrogate._input_labels == ["x1", "x2"]
    assert alm_obj._surrogate._output_labels == ["z1"]
    assert alm_obj._surrogate._input_bounds == {
        "x1": (-1.5, 1.5), "x2": (-1.5, 1.5)}

    # Check populating a block to finish workflow
    blk = SurrogateBlock(concrete=True)

    blk.build_model(alm_obj._surrogate)

    assert isinstance(blk._inputs, Var)
    assert blk._inputs["x1"].bounds == (-1.5, 1.5)
    assert blk._inputs["x2"].bounds == (-1.5, 1.5)
    assert isinstance(blk._outputs, Var)
    assert blk._outputs["z1"].bounds == (None, None)
    assert isinstance(blk.alamo_constraint, Constraint)
    assert len(blk.alamo_constraint) == 1
