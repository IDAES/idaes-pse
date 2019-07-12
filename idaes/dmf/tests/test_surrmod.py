##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Test surrogate modeling.
"""
# stdlib
import logging
import sys

# third-party
import pytest

# local
from idaes.dmf import surrmod, errors
from idaes.dmf.experiment import Experiment

# for testing
from .util import init_logging, tmp_dmf

__author__ = "Dan Gunter <dkgunter@lbl.gov>"

if sys.platform.startswith("win"):
    pytest.skip("skipping DMF tests on Windows", allow_module_level=True)

init_logging()
_log = logging.getLogger(__name__)


@pytest.fixture
def model_data():
    return {"scoops": [1, 2, 3, 4], "cost": [0.5, 0.9, 1.2, 1.4]}


# Real ALAMO


@pytest.mark.skipif(not surrmod.alamo_enabled, reason="ALAMO is disabled")
def test_init(tmp_dmf):
    _test_init(tmp_dmf)


def _test_init(tmp_dmf):
    _ = surrmod.SurrogateModel(Experiment(tmp_dmf))


@pytest.mark.skipif(not surrmod.alamo_enabled, reason="ALAMO is disabled")
def test_run(tmp_dmf, model_data):
    _test_run(tmp_dmf, model_data)


def _test_run(tmp_dmf, model_data):
    sml = surrmod.SurrogateModel(Experiment(tmp_dmf))
    sml.set_input_data(model_data, ["scoops"], "cost")
    try:
        _ = sml.run()
    except errors.AlamoError:
        raise


# Mock ALAMO


class AlamoMock(object):
    def __init__(self):
        self.doalamo_called = False

    def doalamo(self, x, z, *args, **kwargs):
        self.doalamo_called = True
        self.xdata = x
        self.zdata = z
        return {}

    def __enter__(self):
        """Monkey-patch 'surrmod' to use this class
        """
        surrmod.alamo_enabled = True
        surrmod.alamopy = self  # redirect calls to alamopy.doalamo()
        return self

    def __exit__(self, *args):
        surrmod.alamo_enabled = False
        surrmod.alamopy = None


@pytest.mark.skipif(surrmod.alamo_enabled, reason="ALAMO works, no mocking")
def test_init_mock(tmp_dmf):
    with AlamoMock() as mock:
        _test_init(tmp_dmf)
        assert not mock.doalamo_called


@pytest.mark.skipif(surrmod.alamo_enabled, reason="ALAMO works, no mocking")
def test_run_mock(tmp_dmf, model_data):
    with AlamoMock() as mock:
        _test_run(tmp_dmf, model_data)
        assert mock.doalamo_called
        assert tuple(mock.xdata) == tuple(model_data["scoops"])
        assert tuple(mock.zdata) == tuple(model_data["cost"])


@pytest.mark.skipif(surrmod.alamo_enabled, reason="ALAMO works, no mocking")
def test_bad_columns_mock(tmp_dmf, model_data):
    with AlamoMock() as mock:
        m = surrmod.SurrogateModel(Experiment(tmp_dmf))
        with pytest.raises(KeyError):
            m.set_input_data(model_data, ["snork"], "cost")
