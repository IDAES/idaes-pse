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

from idaes.surrogate.alamopy_new import Alamopy


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

        # alm_obj._vdata_in = np.array([[2.5], [6.5]])
        # alm_obj._vdata_out = np.array([25])

        return alm_obj

    @pytest.mark.unit
    def test_default(self, alm_obj):
        stream = io.StringIO()

        alm_obj._write_alamo_input_to_stream(stream=stream)

        print(stream.getvalue())
        assert False
