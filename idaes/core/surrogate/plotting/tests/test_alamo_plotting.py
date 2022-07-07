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
Tests for Surrogate Plotting Methods
"""

import pytest
import os
import pandas as pd

from pyomo.common.fileutils import this_file_dir
from pyomo.common.tempfiles import TempfileManager

from idaes.core.surrogate.alamopy import AlamoSurrogate
from idaes.core.surrogate.sampling.data_utils import split_training_validation
from idaes.core.surrogate.plotting.sm_plotter import (
    surrogate_scatter2D,
    surrogate_scatter3D,
    surrogate_parity,
    surrogate_residual,
)


@pytest.fixture
def alamo_surrogate():
    # import surrogates from external JSON

    alamo_surrogate = AlamoSurrogate.load_from_file(
        os.path.join(this_file_dir(), "alamo_surrogate.json")
    )

    return alamo_surrogate


@pytest.fixture
def data_validation():
    # import validation data and create dataframe here as well
    # random data subset will change every time, but should be fine

    csv_data = pd.read_csv(os.path.join(this_file_dir(), "reformer-data.csv"))
    data = csv_data.sample(n=100)  # randomly sample points for validation

    input_data = data.iloc[:, :2]
    input_labels = input_data.columns

    n_data = data[input_labels[0]].size  # size = 100
    data_training, data_validation = split_training_validation(data, 0.8, seed=n_data)

    return data_validation


# --------------------------------------------------------------------------


@pytest.mark.unit
def test_scatter2D_alamo(alamo_surrogate, data_validation):

    with TempfileManager.new_context() as tf:
        # note - a failure 'The process cannot access the file because it is
        # being used by another process' could occur if an internal error
        # arises before the results file is closed inside the surrogate method

        # create and step into new temporary directory
        dname = tf.mkdtemp()
        filename = os.path.join(dname, "results.pdf")
        surrogate_scatter2D(
            alamo_surrogate, data_validation, filename=filename, show=False
        )

        assert os.path.exists(filename)  # PDF results file


@pytest.mark.unit
def test_scatter3D_alamo(alamo_surrogate, data_validation):

    with TempfileManager.new_context() as tf:
        # note - a failure 'The process cannot access the file because it is
        # being used by another process' could occur if an internal error
        # arises before the results file is closed inside the surrogate method

        # create and step into new temporary directory
        dname = tf.mkdtemp()
        filename = os.path.join(dname, "results.pdf")
        surrogate_scatter3D(
            alamo_surrogate, data_validation, filename=filename, show=False
        )

        assert os.path.exists(filename)  # PDF results file


@pytest.mark.unit
def test_parity_alamo(alamo_surrogate, data_validation):

    with TempfileManager.new_context() as tf:
        # note - a failure 'The process cannot access the file because it is
        # being used by another process' could occur if an internal error
        # arises before the results file is closed inside the surrogate method

        # create and step into new temporary directory
        dname = tf.mkdtemp()
        filename = os.path.join(dname, "results.pdf")
        surrogate_parity(
            alamo_surrogate, data_validation, filename=filename, show=False
        )

        assert os.path.exists(filename)  # PDF results file


@pytest.mark.unit
def test_residual_alamo(alamo_surrogate, data_validation):

    with TempfileManager.new_context() as tf:
        # note - a failure 'The process cannot access the file because it is
        # being used by another process' could occur if an internal error
        # arises before the results file is closed inside the surrogate method

        # create and step into new temporary directory
        dname = tf.mkdtemp()
        filename = os.path.join(dname, "results.pdf")
        surrogate_residual(
            alamo_surrogate, data_validation, filename=filename, show=False
        )

        assert os.path.exists(filename)  # PDF results file


@pytest.mark.unit
def test_scatter2D_noPDF_alamo(alamo_surrogate, data_validation):

    with TempfileManager.new_context() as tf:
        # note - a failure 'The process cannot access the file because it is
        # being used by another process' could occur if an internal error
        # arises before the results file is closed inside the surrogate method

        # create and step into new temporary directory
        dname = tf.mkdtemp()
        surrogate_scatter2D(alamo_surrogate, data_validation, show=False)

        for file in list(os.walk(dname)):  # check entire temp directory
            assert file[-4:] != ".pdf"  # no PDF files should be created


@pytest.mark.unit
def test_scatter3D_noPDF_alamo(alamo_surrogate, data_validation):

    with TempfileManager.new_context() as tf:
        # note - a failure 'The process cannot access the file because it is
        # being used by another process' could occur if an internal error
        # arises before the results file is closed inside the surrogate method

        # create and step into new temporary directory
        dname = tf.mkdtemp()
        surrogate_scatter3D(alamo_surrogate, data_validation, show=False)

        for file in list(os.walk(dname)):  # check entire temp directory
            assert file[-4:] != ".pdf"  # no PDF files should be created


@pytest.mark.unit
def test_parity_noPDF_alamo(alamo_surrogate, data_validation):

    with TempfileManager.new_context() as tf:
        # note - a failure 'The process cannot access the file because it is
        # being used by another process' could occur if an internal error
        # arises before the results file is closed inside the surrogate method

        # create and step into new temporary directory
        dname = tf.mkdtemp()
        surrogate_parity(alamo_surrogate, data_validation, show=False)

        for file in list(os.walk(dname)):  # check entire temp directory
            assert file[-4:] != ".pdf"  # no PDF files should be created


@pytest.mark.unit
def test_residual_noPDF_alamo(alamo_surrogate, data_validation):

    with TempfileManager.new_context() as tf:
        # note - a failure 'The process cannot access the file because it is
        # being used by another process' could occur if an internal error
        # arises before the results file is closed inside the surrogate method

        # create and step into new temporary directory
        dname = tf.mkdtemp()
        surrogate_residual(alamo_surrogate, data_validation, show=False)

        for file in list(os.walk(dname)):  # check entire temp directory
            assert file[-4:] != ".pdf"  # no PDF files should be created
