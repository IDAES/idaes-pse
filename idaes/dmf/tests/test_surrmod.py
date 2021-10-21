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
Test surrogate modeling.
"""
# stdlib
import logging
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Union


# third-party
import pytest

# local
from idaes.dmf import surrmod, DMF
from idaes.dmf.experiment import Experiment

# for testing
from .util import init_logging

__author__ = "Dan Gunter"

init_logging()
_log = logging.getLogger(__name__)
scratch_dir: Union[str, None] = None
scratch_path: Union[Path, None] = None


def setup_module(module):
    global scratch_dir, scratch_path
    scratch_dir = TemporaryDirectory(prefix="idaes_dmf_")  # easier to remove later
    scratch_path = Path(scratch_dir.name)


def teardown_module(module):
    global scratch_dir
    del scratch_dir


model_data = {"scoops": [1, 2, 3, 4], "bowl": [2, 4, 6, 7], "cost": [0.5, 0.9, 1.2, 1.4]}


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


@pytest.mark.skipif(not surrmod.alamo_enabled, reason="ALAMO is disabled")
@pytest.mark.unit
def test_init():
    tmp_dir = scratch_path / "init"
    dmf = DMF(path=tmp_dir, create=True)
    surrmod.SurrogateModel(Experiment(dmf))


@pytest.mark.skipif(not surrmod.alamo_enabled, reason="ALAMO is disabled")
@pytest.mark.unit
def test_run():
    tmp_dir = scratch_path / "run"
    dmf = DMF(path=tmp_dir, create=True)
    sml = surrmod.SurrogateModel(Experiment(dmf))
    sml.set_input_data(model_data, ["scoops", "bowl"], "cost")
    sml.run()


@pytest.mark.unit
def test_mock_alamo():
    tmp_dir = scratch_path / "mock_alamo"
    dmf = DMF(path=tmp_dir, create=True)
    with AlamoMock() as mock:
        sml = surrmod.SurrogateModel(Experiment(dmf))
        assert not mock.doalamo_called
        sml.set_input_data(model_data, ["scoops","bowl"], "cost")
        sml.run()
        assert mock.doalamo_called
        # assert tuple(mock.xdata) == tuple(model_data["scoops"],model_data["bowl"])
        assert tuple(mock.zdata) == tuple(model_data["cost"])


@pytest.mark.unit
def test_bad_columns_mock():
    tmp_dir = scratch_path / "bad_columns_mock"
    dmf = DMF(path=tmp_dir, create=True)
    with AlamoMock() as mock:
        m = surrmod.SurrogateModel(Experiment(dmf))
        with pytest.raises(KeyError):
            m.set_input_data(model_data, ["snork"], "cost")
