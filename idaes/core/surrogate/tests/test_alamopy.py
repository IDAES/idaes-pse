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
import pandas as pd
import io
import os
from math import sin, cos, log, exp
from pathlib import Path
from io import StringIO

from pyomo.environ import Var, Constraint
from pyomo.common.tempfiles import TempfileManager

from idaes.core.surrogate.alamopy import (
    AlamoTrainer,
    AlamoSurrogate,
    Modelers,
    Screener,
    alamo,
)
from idaes.core.surrogate.surrogate_block import SurrogateBlock
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.surrogate.metrics import compute_fit_metrics


dirpath = Path(__file__).parent.resolve()


# String representation of json output for testing
jstring = (
    '{"surrogate": {"z1": " z1 == 3.9999999999925446303450 * x1**2 - '
    "4.0000000000020765611453 * x2**2 - "
    "2.0999999999859380039879 * x1**4 + "
    "4.0000000000043112180492 * x2**4 + "
    "0.33333333332782633107172 * x1**6 + "
    '0.99999999999972988273811 * x1*x2"}, '
    '"input_labels": ["x1", "x2"], '
    '"output_labels": ["z1"], '
    '"input_bounds": {"x1": [0, 5], "x2": [0, 10]}}'
)


class TestAlamoTrainer:
    @pytest.fixture
    def alamo_trainer(self):
        data = {"x1": [1, 2, 3, 4], "x2": [5, 6, 7, 8], "z1": [10, 20, 30, 40]}
        data = pd.DataFrame(data)

        input_labels = ["x1", "x2"]
        output_labels = ["z1"]
        bnds = {"x1": (0, 5), "x2": (0, 10)}

        return AlamoTrainer(
            input_labels=input_labels,
            output_labels=output_labels,
            input_bounds=bnds,
            training_dataframe=data,
        )

    @pytest.mark.unit
    def test_get_files(self, alamo_trainer):
        alamo_trainer.config.filename = "foo.alm"
        alamo_trainer._get_files()
        assert alamo_trainer._almfile == "foo.alm"
        assert alamo_trainer._trcfile == "foo.trc"

    @pytest.mark.unit
    def test_get_files_default(self, alamo_trainer):
        alamo_trainer._get_files()

        assert alamo_trainer._almfile is not None
        assert os.path.exists(alamo_trainer._almfile)
        assert (
            str(alamo_trainer._trcfile).split(".")[0]
            == str(alamo_trainer._almfile).split(".")[0]
        )

        alamo_trainer._remove_temp_files()

        # Check that we cleaned-up after ourselves
        assert not os.path.exists(alamo_trainer._almfile)

    @pytest.mark.unit
    def test_get_files_exists(self, alamo_trainer):
        alamo_trainer.config.filename = os.path.join(dirpath, "alamotrace.trc")

        with pytest.raises(
            FileExistsError,
            match="alamotrace.trc already exists. Either choose a new "
            "file name or set overwrite_files = True",
        ):
            alamo_trainer._get_files()

    @pytest.mark.unit
    def test_spaces_in_input(self):
        data = {"x 1": [1, 2, 3, 4], "x2": [5, 6, 7, 8], "z1": [10, 20, 30, 40]}
        data = pd.DataFrame(data)

        input_labels = ["x 1", "x2"]
        output_labels = ["z1"]
        bnds = {"x 1": (0, 5), "x2": (0, 10)}

        alamo_trainer = AlamoTrainer(
            input_labels=input_labels,
            output_labels=output_labels,
            input_bounds=bnds,
            training_dataframe=data,
        )

        with pytest.raises(
            ValueError,
            match="ALAMO does not support the presence of spaces in variable names. "
            "Found space in input names: x 1.",
        ):
            alamo_trainer._verify_inputs()

    @pytest.mark.unit
    def test_spaces_in_output(self):
        data = {"x1": [1, 2, 3, 4], "x2": [5, 6, 7, 8], "z 1": [10, 20, 30, 40]}
        data = pd.DataFrame(data)

        input_labels = ["x1", "x2"]
        output_labels = ["z 1"]
        bnds = {"x1": (0, 5), "x2": (0, 10)}

        alamo_trainer = AlamoTrainer(
            input_labels=input_labels,
            output_labels=output_labels,
            input_bounds=bnds,
            training_dataframe=data,
        )

        with pytest.raises(
            ValueError,
            match="ALAMO does not support the presence of spaces in variable names. "
            "Found space in output names: z 1.",
        ):
            alamo_trainer._verify_inputs()

    @pytest.mark.unit
    def test_writer_default(self, alamo_trainer):
        stream = io.StringIO()
        alamo_trainer._write_alm_to_stream(stream=stream)

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
            "END_DATA\n"
        )

    @pytest.mark.unit
    def test_writer_reordered(self):
        # Test case where training data is out of order and has extra columns
        data = {
            "junk": [100, 200, 300, 400],
            "x1": [1, 2, 3, 4],
            "z1": [10, 20, 30, 40],
            "x2": [5, 6, 7, 8],
        }
        data = pd.DataFrame(data)

        input_labels = ["x1", "x2"]
        output_labels = ["z1"]
        bnds = {"x1": (0, 5), "x2": (0, 10)}

        alamo_trainer = AlamoTrainer(
            input_labels=input_labels,
            output_labels=output_labels,
            input_bounds=bnds,
            training_dataframe=data,
        )

        stream = io.StringIO()
        alamo_trainer._write_alm_to_stream(stream=stream)

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
            "END_DATA\n"
        )

    @pytest.mark.unit
    def test_writer_custom_basis(self, alamo_trainer):
        alamo_trainer.config.custom_basis_functions = ["sin(x1*x2)", "x1*tanh(x2)"]
        stream = io.StringIO()
        alamo_trainer._write_alm_to_stream(stream=stream)

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
            "END_CUSTOMBAS\n"
        )

    @pytest.mark.unit
    def test_writer_min_max_equal(self, alamo_trainer):
        alamo_trainer._input_bounds = {"x1": (0, 5), "x2": (10, 10)}
        stream = io.StringIO()

        with pytest.raises(
            ConfigurationError,
            match="ALAMO configuration error: upper and lower "
            "bounds on input x2 are equal.",
        ):
            alamo_trainer._write_alm_to_stream(stream=stream)

    @pytest.mark.unit
    def test_writer_min_max_reversed(self, alamo_trainer):
        alamo_trainer._input_bounds = {"x1": (0, 5), "x2": (15, 10)}
        stream = io.StringIO()

        with pytest.raises(
            ConfigurationError,
            match="ALAMO configuration error: upper bound is "
            "less than lower bound for input x2.",
        ):
            alamo_trainer._write_alm_to_stream(stream=stream)

    @pytest.mark.unit
    def test_writer_trace(self, alamo_trainer):
        stream = io.StringIO()
        alamo_trainer._write_alm_to_stream(stream=stream, trace_fname="foo.bar")

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
            "END_DATA\n"
        )

    @pytest.mark.unit
    def test_writer_validation_data(self):
        training_data = {"x1": [1, 2, 3, 4], "x2": [5, 6, 7, 8], "z1": [10, 20, 30, 40]}
        training_data = pd.DataFrame(training_data)
        validation_data = {"x1": [2.5], "x2": [6.5], "z1": [25]}
        validation_data = pd.DataFrame(validation_data)
        input_labels = ["x1", "x2"]
        output_labels = ["z1"]
        bnds = {"x1": (0, 5), "x2": (0, 10)}

        alamo_trainer = AlamoTrainer(
            input_labels=input_labels,
            output_labels=output_labels,
            input_bounds=bnds,
            training_dataframe=training_data,
            validation_dataframe=validation_data,
        )

        stream = io.StringIO()
        alamo_trainer._write_alm_to_stream(stream=stream)

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
            "END_VALDATA\n"
        )

    @pytest.mark.unit
    def test_writer_validation_data_reordered(self):
        # Test case where validation data is out of order and has extra columns
        training_data = {"x1": [1, 2, 3, 4], "x2": [5, 6, 7, 8], "z1": [10, 20, 30, 40]}
        training_data = pd.DataFrame(training_data)
        validation_data = {"junk": [100], "z1": [25], "x1": [2.5], "x2": [6.5]}
        validation_data = pd.DataFrame(validation_data)
        input_labels = ["x1", "x2"]
        output_labels = ["z1"]
        bnds = {"x1": (0, 5), "x2": (0, 10)}

        alamo_trainer = AlamoTrainer(
            input_labels=input_labels,
            output_labels=output_labels,
            input_bounds=bnds,
            training_dataframe=training_data,
            validation_dataframe=validation_data,
        )

        stream = io.StringIO()
        alamo_trainer._write_alm_to_stream(stream=stream)

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
            "END_VALDATA\n"
        )

    @pytest.mark.unit
    def test_writer_full_config(self, alamo_trainer):
        alamo_trainer.config.xfactor = [1.2, 1.6]
        alamo_trainer.config.xscaling = True
        alamo_trainer.config.scalez = True

        alamo_trainer.config.monomialpower = [0.5, 2]
        alamo_trainer.config.multi2power = [1, 2, 3]
        alamo_trainer.config.multi3power = [1, 2, 3, 4]
        alamo_trainer.config.ratiopower = [1, 2, 3, 4, 5]

        alamo_trainer.config.constant = True
        alamo_trainer.config.linfcns = False
        alamo_trainer.config.expfcns = True
        alamo_trainer.config.logfcns = False
        alamo_trainer.config.sinfcns = True
        alamo_trainer.config.cosfcns = False
        alamo_trainer.config.grbfcns = True
        alamo_trainer.config.rbfparam = 7

        alamo_trainer.config.modeler = Modelers.AICc
        alamo_trainer.config.builder = True
        alamo_trainer.config.backstepper = False
        alamo_trainer.config.convpen = 42
        alamo_trainer.config.screener = Screener.SIS
        alamo_trainer.config.ncvf = 4
        alamo_trainer.config.sismult = 5

        alamo_trainer.config.maxiter = -50
        alamo_trainer.config.maxtime = 500
        alamo_trainer.config.datalimitterms = True
        alamo_trainer.config.maxterms = [-1, 2]
        alamo_trainer.config.minterms = [2, -1]
        alamo_trainer.config.numlimitbasis = False
        alamo_trainer.config.exclude = [1, 0]
        alamo_trainer.config.ignore = 1
        alamo_trainer.config.xisint = [1, 0]
        alamo_trainer.config.zisint = 0

        alamo_trainer.config.tolrelmetric = 8.1
        alamo_trainer.config.tolabsmetric = 9.2
        alamo_trainer.config.tolmeanerror = 10.3
        alamo_trainer.config.tolsse = 11.4

        alamo_trainer.config.mipoptca = 12.5
        alamo_trainer.config.mipoptcr = 13.6
        alamo_trainer.config.linearerror = False
        alamo_trainer.config.GAMS = "foo"
        alamo_trainer.config.GAMSSOLVER = "bar"
        alamo_trainer.config.solvemip = False
        alamo_trainer.config.print_to_screen = True

        stream = io.StringIO()
        alamo_trainer._write_alm_to_stream(stream=stream)

        assert stream.getvalue() == (
            "# IDAES Alamopy input file\n"
            "NINPUTS 2\n"
            "NOUTPUTS 1\n"
            "XLABELS x1 x2\n"
            "ZLABELS z1\n"
            "XMIN 0 0\n"
            "XMAX 5 10\n"
            "NDATA 4\n\n"
            "MONO 2\n"
            "MULTI2 3\n"
            "MULTI3 4\n"
            "RATIOS 5\n"
            "xfactor 1.2 1.6\n"
            "xscaling 1\n"
            "scalez 1\n"
            "monomialpower 0.5 2\n"
            "multi2power 1.0 2.0 3.0\n"
            "multi3power 1.0 2.0 3.0 4.0\n"
            "ratiopower 1.0 2.0 3.0 4.0 5.0\n"
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
            "END_DATA\n"
        )

    @pytest.mark.component
    def test_file_writer(self, alamo_trainer):
        alamo_trainer._get_files()
        alamo_trainer._write_alm_file()

        with open(alamo_trainer._almfile, "r") as f:
            fcont = f.read()
        f.close()
        alamo_trainer._remove_temp_files()

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
            f"TRACEFNAME {alamo_trainer._trcfile}\n\n"
            f"BEGIN_DATA\n"
            f"1 5 10\n"
            f"2 6 20\n"
            f"3 7 30\n"
            f"4 8 40\n"
            f"END_DATA\n"
        )

        # Check for clean up
        assert not os.path.exists(alamo_trainer._almfile)

    @pytest.mark.unit
    @pytest.mark.skipif(not alamo.available(), reason="ALAMO not available")
    def test_call_alamo_empty(self, alamo_trainer, caplog):
        alamo_trainer._almfile = "test"

        # To be safe, clear the logger record before test
        caplog.clear()

        alamo_trainer._call_alamo()

        for record in caplog.records:
            assert record.levelname == "WARNING"
        assert (
            "ALAMO executable returned non-zero return code. Check the "
            "ALAMO output for more information."
        ) in caplog.text

        assert alamo_trainer._temp_context is None

    @pytest.mark.component
    @pytest.mark.skipif(not alamo.available(), reason="ALAMO not available")
    def test_call_alamo_w_input(self, alamo_trainer):
        cwd = os.getcwd()

        alamo_trainer._temp_context = TempfileManager.new_context()
        alamo_trainer._almfile = os.path.join(dirpath, "alamo_test.alm")
        alamo_trainer._wrkdir = dirpath

        alamo_trainer._temp_context.add_tempfile(
            os.path.join(dirpath, "GUI_trace.trc"), exists=False
        )
        alamo_trainer._temp_context.add_tempfile(
            os.path.join(dirpath, "alamo_test.lst"), exists=False
        )
        rc, almlog = alamo_trainer._call_alamo()
        assert rc == 0
        assert "Normal termination" in almlog

        # Check for clean up
        alamo_trainer._temp_context.release(remove=True)
        assert not os.path.exists(os.path.join(dirpath, "GUI_trace.trc"))
        assert not os.path.exists(os.path.join(dirpath, "alamo_test.lst"))
        assert cwd == os.getcwd()

    @pytest.mark.unit
    def test_read_trace_single(self, alamo_trainer):
        alamo_trainer._trcfile = os.path.join(dirpath, "alamotrace.trc")
        trc = alamo_trainer._read_trace_file(alamo_trainer._trcfile)

        mdict = {
            "z1": (
                " z1 == 3.9999999999925446303450 * x1**2 - "
                "4.0000000000020765611453 * x2**2 - "
                "2.0999999999859380039879 * x1**4 + "
                "4.0000000000043112180492 * x2**4 + "
                "0.33333333332782633107172 * x1**6 + "
                "0.99999999999972988273811 * x1*x2"
            )
        }

        assert trc == {
            "filename": "GUI_camel6.alm",
            "NINPUTS": "2",
            "NOUTPUTS": "1",
            "INITIALPOINTS": "10",
            "OUTPUT": {"z1": "1"},
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
            "SSEOLR": {"z1": "0.602E-28"},
            "SSE": {"z1": "0.976E-23"},
            "RMSE": {"z1": "0.988E-12"},
            "R2": {"z1": "1.00"},
            "ModelSize": {"z1": "6"},
            "BIC": {"z1": "-539."},
            "RIC": {"z1": "32.5"},
            "Cp": {"z1": "2.00"},
            "AICc": {"z1": "-513."},
            "HQC": {"z1": "-543."},
            "MSE": {"z1": "0.325E-23"},
            "SSEp": {"z1": "0.976E-23"},
            "MADp": {"z1": "0.115E-07"},
            "OLRTime": {"z1": "0.46875000E-01"},
            "numOLRs": {"z1": "16384"},
            "OLRoneCalls": {"z1": "15"},
            "OLRoneFails": {"z1": "0"},
            "OLRgsiCalls": {"z1": "0"},
            "OLRgsiFails": {"z1": "0"},
            "OLRdgelCalls": {"z1": "16369"},
            "OLRdgelFails": {"z1": "0"},
            "OLRclrCalls": {"z1": "0"},
            "OLRclrFails": {"z1": "0"},
            "OLRgmsCalls": {"z1": "0"},
            "OLRgmsFails": {"z1": "0"},
            "CLRTime": {"z1": "0.0000000"},
            "numCLRs": {"z1": "0"},
            "MIPTime": {"z1": "0.0000000"},
            "NumMIPs": {"z1": "0"},
            "LassoTime": {"z1": "0.0000000"},
            "Metric1Lasso": {"z1": "0.17976931+309"},
            "Metric2Lasso": {"z1": "0.17976931+309"},
            "LassoSuccess": {"z1": "F"},
            "LassoRed": {"z1": "0.0000000"},
            "nBasInitAct": {"z1": "15"},
            "nBas": {"z1": "15"},
            "SimTime": {"z1": "0.0000000"},
            "SimData": {"z1": "0"},
            "TotData": {"z1": "10"},
            "NdataConv": {"z1": "0"},
            "OtherTime": {"z1": "0.46875000E-01"},
            "NumIters": {"z1": "1"},
            "IterConv": {"z1": "0"},
            "TimeConv": {"z1": "0.0000000"},
            "Step0Time": {"z1": "0.0000000"},
            "Step1Time": {"z1": "0.93750000E-01"},
            "Step2Time": {"z1": "0.0000000"},
            "TotalTime": {"z1": "0.93750000E-01"},
            "AlamoStatus": {"z1": "0"},
            "AlamoVersion": {"z1": "2021.5.8"},
            "Model": mdict,
        }

    @pytest.mark.unit
    def test_read_trace_multi(self, alamo_trainer):
        alamo_trainer._output_labels = ["z1", "z2"]

        alamo_trainer._trcfile = os.path.join(dirpath, "alamotrace2.trc")
        trc = alamo_trainer._read_trace_file(alamo_trainer._trcfile)

        mdict = {
            "z1": (
                " z1 == 3.9999999999925446303450 * x1**2 - "
                "4.0000000000020765611453 * x2**2 - "
                "2.0999999999859380039879 * x1**4 + "
                "4.0000000000043112180492 * x2**4 + "
                "0.33333333332782633107172 * x1**6 + "
                "0.99999999999972988273811 * x1*x2"
            ),
            "z2": (
                " z2 == 0.72267799849202937756409E-001 * x1 + "
                "0.68451684753912708791823E-001 * x2 + "
                "1.0677896915911471165117 * x1**2 - "
                "0.70576464806224348258468 * x2**2 - "
                "0.40286283566554989543640E-001 * x1**3 + "
                "0.67785668021684807038607E-002 * x2**5 - "
                "0.14017881460354553180281 * x1**6 + "
                "0.77207049441576647286212 * x2**6 + "
                "0.42143309951518070910481 * x1*x2 - "
                "0.41818729807213093907503E-001"
            ),
        }

        assert trc == {
            "filename": "GUI_camel6_twooutputs.alm",
            "NINPUTS": "2",
            "NOUTPUTS": "2",
            "INITIALPOINTS": "10",
            "OUTPUT": {"z1": "1", "z2": "2"},
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
            "SSEOLR": {"z1": "0.602E-28", "z2": "0.280E-30"},
            "SSE": {"z1": "0.976E-23", "z2": "0.280E-30"},
            "RMSE": {"z1": "0.988E-12", "z2": "0.167E-15"},
            "R2": {"z1": "1.00", "z2": "1.00"},
            "ModelSize": {"z1": "6", "z2": "10"},
            "BIC": {"z1": "-539.", "z2": "-704."},
            "RIC": {"z1": "32.5", "z2": "54.2"},
            "Cp": {"z1": "2.00", "z2": "10.0"},
            "AICc": {"z1": "-513.", "z2": "-707."},
            "HQC": {"z1": "-543.", "z2": "-710."},
            "MSE": {"z1": "0.325E-23", "z2": "0.280E-30"},
            "SSEp": {"z1": "0.976E-23", "z2": "0.280E-30"},
            "MADp": {"z1": "0.115E-07", "z2": "0.118E-11"},
            "OLRTime": {"z1": "0.15625000E-01", "z2": "0.93750000E-01"},
            "numOLRs": {"z1": "16384", "z2": "30827"},
            "OLRoneCalls": {"z1": "30", "z2": "30"},
            "OLRoneFails": {"z1": "0", "z2": "0"},
            "OLRgsiCalls": {"z1": "0", "z2": "0"},
            "OLRgsiFails": {"z1": "0", "z2": "0"},
            "OLRdgelCalls": {"z1": "47181", "z2": "47181"},
            "OLRdgelFails": {"z1": "0", "z2": "0"},
            "OLRclrCalls": {"z1": "0", "z2": "0"},
            "OLRclrFails": {"z1": "0", "z2": "0"},
            "OLRgmsCalls": {"z1": "0", "z2": "0"},
            "OLRgmsFails": {"z1": "0", "z2": "0"},
            "CLRTime": {"z1": "0.0000000", "z2": "0.0000000"},
            "numCLRs": {"z1": "0", "z2": "0"},
            "MIPTime": {"z1": "0.0000000", "z2": "0.0000000"},
            "NumMIPs": {"z1": "0", "z2": "0"},
            "LassoTime": {"z1": "0.0000000", "z2": "0.0000000"},
            "Metric1Lasso": {"z1": "0.17976931+309", "z2": "0.17976931+309"},
            "Metric2Lasso": {"z1": "0.17976931+309", "z2": "0.17976931+309"},
            "LassoSuccess": {"z1": "F", "z2": "F"},
            "LassoRed": {"z1": "0.0000000", "z2": "0.0000000"},
            "nBasInitAct": {"z1": "15", "z2": "15"},
            "nBas": {"z1": "15", "z2": "15"},
            "SimTime": {"z1": "0.0000000", "z2": "0.0000000"},
            "SimData": {"z1": "0", "z2": "0"},
            "TotData": {"z1": "10", "z2": "10"},
            "NdataConv": {"z1": "0", "z2": "0"},
            "OtherTime": {"z1": "0.31250000E-01", "z2": "0.31250000E-01"},
            "NumIters": {"z1": "1", "z2": "1"},
            "IterConv": {"z1": "0", "z2": "0"},
            "TimeConv": {"z1": "0.0000000", "z2": "0.0000000"},
            "Step0Time": {"z1": "0.0000000", "z2": "0.0000000"},
            "Step1Time": {"z1": "0.15625000", "z2": "0.15625000"},
            "Step2Time": {"z1": "0.0000000", "z2": "0.0000000"},
            "TotalTime": {"z1": "0.46875000E-01", "z2": "0.10937500"},
            "AlamoStatus": {"z1": "0", "z2": "0"},
            "AlamoVersion": {"z1": "2021.5.8", "z2": "2021.5.8"},
            "Model": mdict,
        }

    @pytest.mark.unit
    def test_read_trace_number_mismatch(self, alamo_trainer):
        alamo_trainer._trcfile = os.path.join(dirpath, "alamotrace2.trc")
        with pytest.raises(
            RuntimeError,
            match="Mismatch when reading ALAMO trace file. "
            "Expected OUTPUT = 1, found 2.",
        ):
            alamo_trainer._read_trace_file(alamo_trainer._trcfile)

    @pytest.mark.unit
    def test_read_trace_label_mismatch(self, alamo_trainer):
        alamo_trainer._output_labels = ["z1", "z3"]

        alamo_trainer._trcfile = os.path.join(dirpath, "alamotrace2.trc")
        with pytest.raises(
            RuntimeError,
            match="Mismatch when reading ALAMO trace file. "
            "Label of output variable in expression "
            "\(z2\) does not match expected label \(z3\).",
        ):
            alamo_trainer._read_trace_file(alamo_trainer._trcfile)

    @pytest.mark.unit
    def test_read_trace_w_validation(self, alamo_trainer):
        alamo_trainer._output_labels = ["z1", "z2"]

        alamo_trainer._trcfile = os.path.join(dirpath, "alamotrace_w_validation.trc")
        trc = alamo_trainer._read_trace_file(
            alamo_trainer._trcfile, has_validation_data=True
        )

        mdict = {
            "z1": (
                " z1 ==  - 0.99999999999999689137553 * x1 - "
                "0.20616513343350446289142E-015 * x2 + "
                "1.0000000000000002220446 * x1**2 - x2**2 + "
                "0.99999999999999811262086"
            ),
            "z2": (
                " z2 == 1.0000000000000022204460 * x2 - "
                "1.0000000000000002220446 * x1**2 + x2**2 - "
                "0.72986161565720566182285E-016 * x1**3 - "
                "1.0000000000000091038288"
            ),
        }

        assert trc == {
            "filename": "alamo_error.alm",
            "NINPUTS": "2",
            "NOUTPUTS": "2",
            "INITIALPOINTS": "7",
            "OUTPUT": {"z1": "1", "z2": "2"},
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
            "SSEOLR": {"z1": "0.347E-27", "z2": "0.247E-27"},
            "SSE": {"z1": "0.347E-27", "z2": "0.247E-27"},
            "RMSE": {"z1": "0.704E-14", "z2": "0.594E-14"},
            "R2": {"z1": "1.00", "z2": "1.00"},
            "ModelSize": {"z1": "5", "z2": "5"},
            "BIC": {"z1": "-446.", "z2": "-449."},
            "RIC": {"z1": "20.8", "z2": "20.8"},
            "Cp": {"z1": "3.00", "z2": "3.00"},
            "AICc": {"z1": "-386.", "z2": "-389."},
            "HQC": {"z1": "-450.", "z2": "-452."},
            "MSE": {"z1": "0.347E-27", "z2": "0.247E-27"},
            "SSEp": {"z1": "0.347E-27", "z2": "0.247E-27"},
            "MADp": {"z1": "0.888E-13", "z2": "0.888E-13"},
            "OLRTime": {"z1": "0.0000000", "z2": "0.0000000"},
            "numOLRs": {"z1": "247", "z2": "247"},
            "OLRoneCalls": {"z1": "16", "z2": "16"},
            "OLRoneFails": {"z1": "0", "z2": "0"},
            "OLRgsiCalls": {"z1": "0", "z2": "0"},
            "OLRgsiFails": {"z1": "0", "z2": "0"},
            "OLRdgelCalls": {"z1": "478", "z2": "478"},
            "OLRdgelFails": {"z1": "0", "z2": "0"},
            "OLRclrCalls": {"z1": "0", "z2": "0"},
            "OLRclrFails": {"z1": "0", "z2": "0"},
            "OLRgmsCalls": {"z1": "0", "z2": "0"},
            "OLRgmsFails": {"z1": "0", "z2": "0"},
            "CLRTime": {"z1": "0.0000000", "z2": "0.0000000"},
            "numCLRs": {"z1": "0", "z2": "0"},
            "MIPTime": {"z1": "0.0000000", "z2": "0.0000000"},
            "NumMIPs": {"z1": "0", "z2": "0"},
            "LassoTime": {"z1": "0.0000000", "z2": "0.0000000"},
            "Metric1Lasso": {"z1": "0.17976931+309", "z2": "0.17976931+309"},
            "Metric2Lasso": {"z1": "0.17976931+309", "z2": "0.17976931+309"},
            "LassoSuccess": {"z1": "F", "z2": "F"},
            "LassoRed": {"z1": "0.0000000", "z2": "0.0000000"},
            "nBasInitAct": {"z1": "8", "z2": "8"},
            "nBas": {"z1": "8", "z2": "8"},
            "SimTime": {"z1": "0.0000000", "z2": "0.0000000"},
            "SimData": {"z1": "0", "z2": "0"},
            "TotData": {"z1": "7", "z2": "7"},
            "NdataConv": {"z1": "0", "z2": "0"},
            "OtherTime": {"z1": "0.18750000", "z2": "0.18750000"},
            "NumIters": {"z1": "1", "z2": "1"},
            "IterConv": {"z1": "0", "z2": "0"},
            "TimeConv": {"z1": "0.0000000", "z2": "0.0000000"},
            "Step0Time": {"z1": "0.0000000", "z2": "0.0000000"},
            "Step1Time": {"z1": "0.15625000E-01", "z2": "0.15625000E-01"},
            "Step2Time": {"z1": "0.0000000", "z2": "0.0000000"},
            "TotalTime": {"z1": "0.0000000", "z2": "0.15625000E-01"},
            "AlamoStatus": {"z1": "0", "z2": "0"},
            "AlamoVersion": {"z1": "2021.12.28", "z2": "2021.12.28"},
            "Model": mdict,
        }

    @pytest.mark.unit
    def test_read_trace_invalid_file(self, alamo_trainer):
        alamo_trainer._output_labels = ["z1", "z2"]

        alamo_trainer._trcfile = os.path.join(dirpath, "foo.trc")

        with pytest.raises(
            FileNotFoundError,
            match="Error occured when trying to read the ALAMO trace file - "
            "this probably indicates that a trace file was not created by "
            "the ALAMO executable. Please check the ALAMO output logs.",
        ):
            alamo_trainer._read_trace_file(
                alamo_trainer._trcfile, has_validation_data=True
            )

    @pytest.mark.unit
    def test_populate_results(self, alamo_trainer):
        alamo_trainer._trcfile = os.path.join(dirpath, "alamotrace.trc")
        trc = alamo_trainer._read_trace_file(alamo_trainer._trcfile)
        alamo_trainer._populate_results(trc)

        mdict = {
            "z1": (
                " z1 == 3.9999999999925446303450 * x1**2 - "
                "4.0000000000020765611453 * x2**2 - "
                "2.0999999999859380039879 * x1**4 + "
                "4.0000000000043112180492 * x2**4 + "
                "0.33333333332782633107172 * x1**6 + "
                "0.99999999999972988273811 * x1*x2"
            )
        }

        assert alamo_trainer._results == {
            "filename": "GUI_camel6.alm",
            "NINPUTS": "2",
            "NOUTPUTS": "1",
            "INITIALPOINTS": "10",
            "OUTPUT": {"z1": "1"},
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
            "SSEOLR": {"z1": "0.602E-28"},
            "SSE": {"z1": "0.976E-23"},
            "RMSE": {"z1": "0.988E-12"},
            "R2": {"z1": "1.00"},
            "ModelSize": {"z1": "6"},
            "BIC": {"z1": "-539."},
            "RIC": {"z1": "32.5"},
            "Cp": {"z1": "2.00"},
            "AICc": {"z1": "-513."},
            "HQC": {"z1": "-543."},
            "MSE": {"z1": "0.325E-23"},
            "SSEp": {"z1": "0.976E-23"},
            "MADp": {"z1": "0.115E-07"},
            "OLRTime": {"z1": "0.46875000E-01"},
            "numOLRs": {"z1": "16384"},
            "OLRoneCalls": {"z1": "15"},
            "OLRoneFails": {"z1": "0"},
            "OLRgsiCalls": {"z1": "0"},
            "OLRgsiFails": {"z1": "0"},
            "OLRdgelCalls": {"z1": "16369"},
            "OLRdgelFails": {"z1": "0"},
            "OLRclrCalls": {"z1": "0"},
            "OLRclrFails": {"z1": "0"},
            "OLRgmsCalls": {"z1": "0"},
            "OLRgmsFails": {"z1": "0"},
            "CLRTime": {"z1": "0.0000000"},
            "numCLRs": {"z1": "0"},
            "MIPTime": {"z1": "0.0000000"},
            "NumMIPs": {"z1": "0"},
            "LassoTime": {"z1": "0.0000000"},
            "Metric1Lasso": {"z1": "0.17976931+309"},
            "Metric2Lasso": {"z1": "0.17976931+309"},
            "LassoSuccess": {"z1": "F"},
            "LassoRed": {"z1": "0.0000000"},
            "nBasInitAct": {"z1": "15"},
            "nBas": {"z1": "15"},
            "SimTime": {"z1": "0.0000000"},
            "SimData": {"z1": "0"},
            "TotData": {"z1": "10"},
            "NdataConv": {"z1": "0"},
            "OtherTime": {"z1": "0.46875000E-01"},
            "NumIters": {"z1": "1"},
            "IterConv": {"z1": "0"},
            "TimeConv": {"z1": "0.0000000"},
            "Step0Time": {"z1": "0.0000000"},
            "Step1Time": {"z1": "0.93750000E-01"},
            "Step2Time": {"z1": "0.0000000"},
            "TotalTime": {"z1": "0.93750000E-01"},
            "AlamoStatus": {"z1": "0"},
            "AlamoVersion": {"z1": "2021.5.8"},
            "Model": mdict,
        }

    @pytest.mark.unit
    def test_build_surrogate_object(self, alamo_trainer):
        alamo_trainer._trcfile = os.path.join(dirpath, "alamotrace.trc")
        trc = alamo_trainer._read_trace_file(alamo_trainer._trcfile)
        alamo_trainer._populate_results(trc)
        alamo_object = alamo_trainer._build_surrogate_object()

        assert isinstance(alamo_object, AlamoSurrogate)
        assert alamo_object._surrogate_expressions == {
            "z1": (
                " z1 == 3.9999999999925446303450 * x1**2 - "
                "4.0000000000020765611453 * x2**2 - "
                "2.0999999999859380039879 * x1**4 + "
                "4.0000000000043112180492 * x2**4 + "
                "0.33333333332782633107172 * x1**6 + "
                "0.99999999999972988273811 * x1*x2"
            )
        }
        assert alamo_object._input_labels == ["x1", "x2"]
        assert alamo_object._output_labels == ["z1"]
        assert alamo_object._input_bounds == {"x1": (0, 5), "x2": (0, 10)}


class TestAlamoSurrogate:
    @pytest.fixture
    def alm_surr1(self):
        surrogate_expressions = {
            "z1": (
                " z1 == 3.9999999999925446303450 * x1**2 - "
                "4.0000000000020765611453 * x2**2 - "
                "2.0999999999859380039879 * x1**4 + "
                "4.0000000000043112180492 * x2**4 + "
                "0.33333333332782633107172 * x1**6 + "
                "0.99999999999972988273811 * x1*x2"
            )
        }
        input_labels = ["x1", "x2"]
        output_labels = ["z1"]
        input_bounds = {"x1": (0, 5), "x2": (0, 10)}

        alm_surr1 = AlamoSurrogate(
            surrogate_expressions, input_labels, output_labels, input_bounds
        )

        return alm_surr1

    @pytest.mark.unit
    def test_evaluate_surrogate(self, alm_surr1):
        x = [
            -2,
            -1.8,
            -1.6,
            -1.4,
            -1.2,
            -1.0,
            -0.8,
            -0.6,
            -0.4,
            -0.2,
            0,
            0.2,
            0.4,
            0.6,
            0.8,
            1.0,
            1.2,
            1.4,
            1.6,
            1.8,
            2.0,
        ]

        inputs = np.array([np.tile(x, len(x)), np.repeat(x, len(x))])
        inputs = pd.DataFrame(inputs.transpose(), columns=["x1", "x2"])

        out = alm_surr1.evaluate_surrogate(inputs)
        for i in range(inputs.shape[0]):
            assert pytest.approx(out["z1"][i], rel=1e-8) == (
                3.9999999999925446303450 * inputs["x1"][i] ** 2
                - 4.0000000000020765611453 * inputs["x2"][i] ** 2
                - 2.0999999999859380039879 * inputs["x1"][i] ** 4
                + 4.0000000000043112180492 * inputs["x2"][i] ** 4
                + 0.33333333332782633107172 * inputs["x1"][i] ** 6
                + 0.99999999999972988273811 * inputs["x1"][i] * inputs["x2"][i]
            )

        # Recall evaluate_surrogate to make sure cached fcn works
        out2 = alm_surr1.evaluate_surrogate(inputs)
        for i in range(inputs.shape[0]):
            assert pytest.approx(out2["z1"][i], rel=1e-8) == (
                3.9999999999925446303450 * inputs["x1"][i] ** 2
                - 4.0000000000020765611453 * inputs["x2"][i] ** 2
                - 2.0999999999859380039879 * inputs["x1"][i] ** 4
                + 4.0000000000043112180492 * inputs["x2"][i] ** 4
                + 0.33333333332782633107172 * inputs["x1"][i] ** 6
                + 0.99999999999972988273811 * inputs["x1"][i] * inputs["x2"][i]
            )

    @pytest.mark.unit
    def test_populate_block(self, alm_surr1):
        blk = SurrogateBlock(concrete=True)

        blk.build_model(alm_surr1)
        blk.display()

        assert isinstance(blk.inputs, Var)
        assert blk.inputs["x1"].bounds == (0, 5)
        assert blk.inputs["x2"].bounds == (0, 10)
        assert isinstance(blk.outputs, Var)
        assert blk.outputs["z1"].bounds == (None, None)
        assert isinstance(blk.alamo_constraint, Constraint)
        assert len(blk.alamo_constraint) == 1
        assert str(blk.alamo_constraint["z1"].body) == (
            "outputs[z1] - ("
            "3.9999999999925446*inputs[x1]**2 - "
            "4.000000000002077*inputs[x2]**2 - "
            "2.099999999985938*inputs[x1]**4 + "
            "4.000000000004311*inputs[x2]**4 + "
            "0.33333333332782633*inputs[x1]**6 + "
            "0.9999999999997299*inputs[x1]*inputs[x2])"
        )

    @pytest.fixture
    def alm_surr2(self):
        surrogate_expressions = {
            "z1": (
                "z1 == 3.9999999999925446303450 * x1**2 - "
                "4.0000000000020765611453 * x2**2 - "
                "2.0999999999859380039879 * x1**4 + "
                "4.0000000000043112180492 * x2**4 + "
                "0.33333333332782633107172 * x1**6 + "
                "0.99999999999972988273811 * x1*x2"
            ),
            "z2": (
                "z2 == 0.72267799849202937756409E-001 * x1 + "
                "0.68451684753912708791823E-001 * x2 + "
                "1.0677896915911471165117 * x1**2 - "
                "0.70576464806224348258468 * x2**2 - "
                "0.40286283566554989543640E-001 * x1**3 + "
                "0.67785668021684807038607E-002 * x2**5 - "
                "0.14017881460354553180281 * x1**6 + "
                "0.77207049441576647286212 * x2**6 + "
                "0.42143309951518070910481 * x1*x2 - "
                "0.41818729807213093907503E-001"
            ),
        }
        input_labels = ["x1", "x2"]
        output_labels = ["z1", "z2"]

        alm_surr2 = AlamoSurrogate(surrogate_expressions, input_labels, output_labels)

        return alm_surr2

    @pytest.mark.unit
    def test_evaluate_surrogate_multi(self, alm_surr2):
        x = [
            -2,
            -1.8,
            -1.6,
            -1.4,
            -1.2,
            -1.0,
            -0.8,
            -0.6,
            -0.4,
            -0.2,
            0,
            0.2,
            0.4,
            0.6,
            0.8,
            1.0,
            1.2,
            1.4,
            1.6,
            1.8,
            2.0,
        ]

        inputs = np.array([np.tile(x, len(x)), np.repeat(x, len(x))])
        inputs = pd.DataFrame(inputs.transpose(), columns=["x1", "x2"])

        out = alm_surr2.evaluate_surrogate(inputs)
        for i in range(inputs.shape[0]):
            assert pytest.approx(out["z1"][i], rel=1e-8) == (
                3.9999999999925446303450 * inputs["x1"][i] ** 2
                - 4.0000000000020765611453 * inputs["x2"][i] ** 2
                - 2.0999999999859380039879 * inputs["x1"][i] ** 4
                + 4.0000000000043112180492 * inputs["x2"][i] ** 4
                + 0.33333333332782633107172 * inputs["x1"][i] ** 6
                + 0.99999999999972988273811 * inputs["x1"][i] * inputs["x2"][i]
            )
            assert pytest.approx(out["z2"][i], rel=1e-8) == (
                0.72267799849202937756409e-001 * inputs["x1"][i]
                + 0.68451684753912708791823e-001 * inputs["x2"][i]
                + 1.0677896915911471165117 * inputs["x1"][i] ** 2
                - 0.70576464806224348258468 * inputs["x2"][i] ** 2
                - 0.40286283566554989543640e-001 * inputs["x1"][i] ** 3
                + 0.67785668021684807038607e-002 * inputs["x2"][i] ** 5
                - 0.14017881460354553180281 * inputs["x1"][i] ** 6
                + 0.77207049441576647286212 * inputs["x2"][i] ** 6
                + 0.42143309951518070910481 * inputs["x1"][i] * inputs["x2"][i]
                - 0.41818729807213093907503e-001
            )

        # Recall evaluate_surrogate to make sure cached fcn works
        out2 = alm_surr2.evaluate_surrogate(inputs)
        for i in range(inputs.shape[0]):
            assert pytest.approx(out2["z1"][i], rel=1e-8) == (
                3.9999999999925446303450 * inputs["x1"][i] ** 2
                - 4.0000000000020765611453 * inputs["x2"][i] ** 2
                - 2.0999999999859380039879 * inputs["x1"][i] ** 4
                + 4.0000000000043112180492 * inputs["x2"][i] ** 4
                + 0.33333333332782633107172 * inputs["x1"][i] ** 6
                + 0.99999999999972988273811 * inputs["x1"][i] * inputs["x2"][i]
            )
            assert pytest.approx(out2["z2"][i], rel=1e-8) == (
                0.72267799849202937756409e-001 * inputs["x1"][i]
                + 0.68451684753912708791823e-001 * inputs["x2"][i]
                + 1.0677896915911471165117 * inputs["x1"][i] ** 2
                - 0.70576464806224348258468 * inputs["x2"][i] ** 2
                - 0.40286283566554989543640e-001 * inputs["x1"][i] ** 3
                + 0.67785668021684807038607e-002 * inputs["x2"][i] ** 5
                - 0.14017881460354553180281 * inputs["x1"][i] ** 6
                + 0.77207049441576647286212 * inputs["x2"][i] ** 6
                + 0.42143309951518070910481 * inputs["x1"][i] * inputs["x2"][i]
                - 0.41818729807213093907503e-001
            )

    @pytest.mark.unit
    def test_populate_block_multi(self, alm_surr2):
        blk = SurrogateBlock(concrete=True)

        blk.build_model(alm_surr2)

        assert isinstance(blk.inputs, Var)
        assert blk.inputs["x1"].bounds == (None, None)
        assert blk.inputs["x2"].bounds == (None, None)
        assert isinstance(blk.outputs, Var)
        assert blk.outputs["z1"].bounds == (None, None)
        assert blk.outputs["z2"].bounds == (None, None)
        assert isinstance(blk.alamo_constraint, Constraint)
        assert len(blk.alamo_constraint) == 2
        assert str(blk.alamo_constraint["z1"].body) == (
            "outputs[z1] - ("
            "3.9999999999925446*inputs[x1]**2 - "
            "4.000000000002077*inputs[x2]**2 - "
            "2.099999999985938*inputs[x1]**4 + "
            "4.000000000004311*inputs[x2]**4 + "
            "0.33333333332782633*inputs[x1]**6 + "
            "0.9999999999997299*inputs[x1]*inputs[x2])"
        )
        assert str(blk.alamo_constraint["z2"].body) == (
            "outputs[z2] - ("
            "0.07226779984920294*inputs[x1] + "
            "0.06845168475391271*inputs[x2] + "
            "1.0677896915911471*inputs[x1]**2 - "
            "0.7057646480622435*inputs[x2]**2 - "
            "0.04028628356655499*inputs[x1]**3 + "
            "0.006778566802168481*inputs[x2]**5 - "
            "0.14017881460354553*inputs[x1]**6 + "
            "0.7720704944157665*inputs[x2]**6 + "
            "0.4214330995151807*inputs[x1]*inputs[x2] - "
            "0.041818729807213094)"
        )

    @pytest.fixture
    def alm_surr3(self):
        surrogate_expressions = {
            "z1": (
                " z1 == 2*sin(x1**2) - 3*cos(x2**3) - " "4*log(x1**4) + 5*exp(x2**5)"
            )
        }
        input_labels = ["x1", "x2"]
        output_labels = ["z1"]

        alm_surr3 = AlamoSurrogate(surrogate_expressions, input_labels, output_labels)

        return alm_surr3

    @pytest.mark.unit
    def test_evaluate_surrogate_funcs(self, alm_surr3):
        x = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]

        inputs = np.array([np.tile(x, len(x)), np.repeat(x, len(x))])
        inputs = pd.DataFrame(inputs.transpose(), columns=["x1", "x2"])

        out = alm_surr3.evaluate_surrogate(inputs)
        for i in range(inputs.shape[0]):
            assert pytest.approx(out["z1"][i], rel=1e-8) == (
                2 * sin(inputs["x1"][i] ** 2)
                - 3 * cos(inputs["x2"][i] ** 3)
                - 4 * log(inputs["x1"][i] ** 4)
                + 5 * exp(inputs["x2"][i] ** 5)
            )

        # Recall evaluate_surrogate to make sure cached fcn works
        out2 = alm_surr3.evaluate_surrogate(inputs)
        for i in range(inputs.shape[0]):
            assert pytest.approx(out2["z1"][i], rel=1e-8) == (
                2 * sin(inputs["x1"][i] ** 2)
                - 3 * cos(inputs["x2"][i] ** 3)
                - 4 * log(inputs["x1"][i] ** 4)
                + 5 * exp(inputs["x2"][i] ** 5)
            )

    @pytest.mark.unit
    def test_populate_block_funcs(self, alm_surr3):
        blk = SurrogateBlock(concrete=True)

        blk.build_model(alm_surr3)

        assert isinstance(blk.inputs, Var)
        assert isinstance(blk.outputs, Var)
        assert isinstance(blk.alamo_constraint, Constraint)
        assert len(blk.alamo_constraint) == 1
        assert str(blk.alamo_constraint["z1"].body) == (
            "outputs[z1] - (2*sin(inputs[x1]**2) - 3*cos(inputs[x2]**3) - "
            "4*log(inputs[x1]**4) + 5*exp(inputs[x2]**5))"
        )

    @pytest.mark.unit
    def test_save(self, alm_surr1):
        stream = StringIO()
        alm_surr1.save(stream)
        assert stream.getvalue() == jstring

    @pytest.mark.unit
    def test_load(self):
        alm_surr = AlamoSurrogate.load(StringIO(jstring))

        assert alm_surr._surrogate_expressions == {
            "z1": " z1 == 3.9999999999925446303450 * x1**2 - "
            "4.0000000000020765611453 * x2**2 - "
            "2.0999999999859380039879 * x1**4 + "
            "4.0000000000043112180492 * x2**4 + "
            "0.33333333332782633107172 * x1**6 + "
            "0.99999999999972988273811 * x1*x2"
        }
        assert alm_surr._input_labels == ["x1", "x2"]
        assert alm_surr._output_labels == ["z1"]
        assert alm_surr._input_bounds == {"x1": (0, 5), "x2": (0, 10)}

    @pytest.mark.unit
    def test_save_load(self, alm_surr1):
        with TempfileManager as tf:
            fname = tf.create_tempfile(suffix=".json")
            alm_surr1.save_to_file(fname, overwrite=True)

            assert os.path.isfile(fname)

            with open(fname, "r") as f:
                js = f.read()
            f.close()

            alm_load = AlamoSurrogate.load_from_file(fname)

        # Check file contents
        assert js == jstring

        # Check loaded object
        assert isinstance(alm_load, AlamoSurrogate)
        assert alm_load._surrogate_expressions == {
            "z1": " z1 == 3.9999999999925446303450 * x1**2 - "
            "4.0000000000020765611453 * x2**2 - "
            "2.0999999999859380039879 * x1**4 + "
            "4.0000000000043112180492 * x2**4 + "
            "0.33333333332782633107172 * x1**6 + "
            "0.99999999999972988273811 * x1*x2"
        }
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
                alm_surr1.save_to_file(fname)

        # Check for clean up
        assert not os.path.isfile(fname)


@pytest.mark.skipif(alamo.executable is None, reason="ALAMO not available")
@pytest.mark.integration
class TestWorkflow:
    training_data = pd.DataFrame(
        np.array(
            [
                [0.353837234435, 0.99275270941666, 0.762878272854],
                [0.904978848612, -0.746908518721, 0.387963718723],
                [0.643706630938, -0.617496599522, -0.0205375902284],
                [1.29881420688, 0.305594881575, 2.43011137696],
                [1.35791650867, 0.351045058258, 2.36989368612],
                [0.938369314089, -0.525167416293, 0.829756159423],
                [-1.46593541641, 0.383902178482, 1.14054797964],
                [-0.374378293218, -0.689730440659, -0.219122783909],
                [0.690326213554, 0.569364994374, 0.982068847698],
                [-0.961163301329, 0.499471920546, 0.936855365038],
            ]
        ),
        columns=["x1", "x2", "z1"],
    )

    @pytest.fixture(scope="class")
    def alamo_trainer(self):
        # Test end-to-end workflow with a simple problem.
        bnds = {"x1": (-1.5, 1.5), "x2": (-1.5, 1.5)}
        alamo_trainer = AlamoTrainer(
            input_labels=["x1", "x2"],
            output_labels=["z1"],
            input_bounds=bnds,
            training_dataframe=TestWorkflow.training_data,
        )

        alamo_trainer.config.linfcns = True
        alamo_trainer.config.monomialpower = [2, 3, 4, 5, 6]
        alamo_trainer.config.multi2power = [1, 2]

        (
            alamo_trainer.status,
            alamo_trainer.alamo_object,
            alamo_trainer.msg,
        ) = alamo_trainer.train_surrogate()

        return alamo_trainer

    def test_execution(self, alamo_trainer):
        # Check execution
        assert alamo_trainer.status is True
        assert alamo_trainer.msg == " Normal termination"

        # Check temp file clean up
        assert alamo_trainer._temp_context is None
        assert not os.path.exists(alamo_trainer._almfile)
        assert not os.path.exists(alamo_trainer._trcfile)

    def test_alamo_results(self, alamo_trainer):
        assert alamo_trainer._results is not None
        assert alamo_trainer._results["NINPUTS"] == "2"
        assert alamo_trainer._results["NOUTPUTS"] == "1"
        assert alamo_trainer._results["SSEOLR"] == {"z1": "0.373E-29"}
        assert alamo_trainer._results["SSE"] == {"z1": "0.976E-23"}
        assert alamo_trainer._results["RMSE"] == {"z1": "0.988E-12"}
        assert alamo_trainer._results["R2"] == {"z1": "1.00"}
        assert alamo_trainer._results["ModelSize"] == {"z1": "6"}
        assert alamo_trainer._results["BIC"] == {"z1": "-539."}
        assert alamo_trainer._results["RIC"] == {"z1": "32.5"}
        assert alamo_trainer._results["Cp"] == {"z1": "2.00"}
        assert alamo_trainer._results["AICc"] == {"z1": "-513."}
        assert alamo_trainer._results["HQC"] == {"z1": "-543."}
        assert alamo_trainer._results["MSE"] == {"z1": "0.325E-23"}
        assert alamo_trainer._results["SSEp"] == {"z1": "0.976E-23"}
        assert alamo_trainer._results["MADp"] == {"z1": "0.115E-07"}

        assert alamo_trainer._results["Model"] == {
            "z1": " z1 == 3.9999999999925432980774 * x1**2 - "
            "4.0000000000020792256805 * x2**2 - "
            "2.0999999999859380039879 * x1**4 + "
            "4.0000000000043085535140 * x2**4 + "
            "0.33333333332782683067208 * x1**6 + "
            "0.99999999999973088193883 * x1*x2"
        }

    def test_alamo_object(self, alamo_trainer):
        alamo_object = alamo_trainer.alamo_object
        assert isinstance(alamo_object, AlamoSurrogate)
        assert alamo_object._surrogate_expressions == {
            "z1": " z1 == 3.9999999999925432980774 * x1**2 - "
            "4.0000000000020792256805 * x2**2 - "
            "2.0999999999859380039879 * x1**4 + "
            "4.0000000000043085535140 * x2**4 + "
            "0.33333333332782683067208 * x1**6 + "
            "0.99999999999973088193883 * x1*x2"
        }
        assert alamo_object._input_labels == ["x1", "x2"]
        assert alamo_object._output_labels == ["z1"]
        assert alamo_object._input_bounds == {"x1": (-1.5, 1.5), "x2": (-1.5, 1.5)}

        # Check populating a block to finish workflow
        blk = SurrogateBlock(concrete=True)

        blk.build_model(alamo_object)

        assert isinstance(blk.inputs, Var)
        assert blk.inputs["x1"].bounds == (-1.5, 1.5)
        assert blk.inputs["x2"].bounds == (-1.5, 1.5)
        assert isinstance(blk.outputs, Var)
        assert blk.outputs["z1"].bounds == (None, None)
        assert isinstance(blk.alamo_constraint, Constraint)
        assert len(blk.alamo_constraint) == 1

    def test_metrics(self, alamo_trainer):
        alamo_object = alamo_trainer.alamo_object

        metrics = compute_fit_metrics(alamo_object, TestWorkflow.training_data)

        assert isinstance(metrics, dict)

        # ALAMO metrics are consistently 10x larger than they should be
        assert metrics["z1"]["SSE"] == pytest.approx(
            float(alamo_trainer._results["SSE"]["z1"]) * 0.1, rel=1e-8
        )
        assert metrics["z1"]["R2"] == pytest.approx(
            float(alamo_trainer._results["R2"]["z1"]), rel=1e-8
        )
        assert metrics["z1"]["MSE"] == pytest.approx(
            float(alamo_trainer._results["RMSE"]["z1"]) ** 2 * 0.1, rel=1e-8
        )
        assert metrics["z1"]["RMSE"] == pytest.approx(
            float(alamo_trainer._results["RMSE"]["z1"]) * 0.1, rel=1e-8
        )


@pytest.mark.skipif(alamo.executable is None, reason="ALAMO not available")
@pytest.mark.integration
class TestWorkflowValidation:
    training_data = pd.DataFrame(
        np.array(
            [
                [0.353837234435, 0.99275270941666, 0.762878272854],
                [0.904978848612, -0.746908518721, 0.387963718723],
                [0.643706630938, -0.617496599522, -0.0205375902284],
                [1.29881420688, 0.305594881575, 2.43011137696],
                [1.35791650867, 0.351045058258, 2.36989368612],
                [0.938369314089, -0.525167416293, 0.829756159423],
                [-1.46593541641, 0.383902178482, 1.14054797964],
                [-0.374378293218, -0.689730440659, -0.219122783909],
                [0.690326213554, 0.569364994374, 0.982068847698],
                [-0.961163301329, 0.499471920546, 0.936855365038],
            ]
        ),
        columns=["x1", "x2", "z1"],
    )

    @pytest.fixture(scope="class")
    def alamo_trainer(self):
        # Test end-to-end workflow with a simple problem.
        # Here we will send the training dataset as validation as well,
        # to confirm that the mechanics work
        bnds = {"x1": (-1.5, 1.5), "x2": (-1.5, 1.5)}
        alamo_trainer = AlamoTrainer(
            input_labels=["x1", "x2"],
            output_labels=["z1"],
            input_bounds=bnds,
            training_dataframe=TestWorkflow.training_data,
            validation_dataframe=TestWorkflow.training_data,
        )

        alamo_trainer.config.linfcns = True
        alamo_trainer.config.monomialpower = [2, 3, 4, 5, 6]
        alamo_trainer.config.multi2power = [1, 2]

        (
            alamo_trainer.status,
            alamo_trainer.alamo_object,
            alamo_trainer.msg,
        ) = alamo_trainer.train_surrogate()

        return alamo_trainer

    def test_execution(self, alamo_trainer):
        # Check execution
        assert alamo_trainer.status is True
        assert alamo_trainer.msg == " Normal termination"

        # Check temp file clean up
        assert alamo_trainer._temp_context is None
        assert not os.path.exists(alamo_trainer._almfile)
        assert not os.path.exists(alamo_trainer._trcfile)

    def test_alamo_results(self, alamo_trainer):
        assert alamo_trainer._results is not None
        assert alamo_trainer._results["NINPUTS"] == "2"
        assert alamo_trainer._results["NOUTPUTS"] == "1"
        assert alamo_trainer._results["SSEOLR"] == {"z1": "0.373E-29"}
        assert alamo_trainer._results["SSE"] == {"z1": "0.976E-23"}
        assert alamo_trainer._results["RMSE"] == {"z1": "0.988E-12"}
        assert alamo_trainer._results["R2"] == {"z1": "1.00"}
        assert alamo_trainer._results["ModelSize"] == {"z1": "6"}
        assert alamo_trainer._results["BIC"] == {"z1": "-539."}
        assert alamo_trainer._results["RIC"] == {"z1": "32.5"}
        assert alamo_trainer._results["Cp"] == {"z1": "2.00"}
        assert alamo_trainer._results["AICc"] == {"z1": "-513."}
        assert alamo_trainer._results["HQC"] == {"z1": "-543."}
        assert alamo_trainer._results["MSE"] == {"z1": "0.325E-23"}
        assert alamo_trainer._results["SSEp"] == {"z1": "0.976E-23"}
        assert alamo_trainer._results["MADp"] == {"z1": "0.115E-07"}

        assert alamo_trainer._results["Model"] == {
            "z1": " z1 == 3.9999999999925432980774 * x1**2 - "
            "4.0000000000020792256805 * x2**2 - "
            "2.0999999999859380039879 * x1**4 + "
            "4.0000000000043085535140 * x2**4 + "
            "0.33333333332782683067208 * x1**6 + "
            "0.99999999999973088193883 * x1*x2"
        }
