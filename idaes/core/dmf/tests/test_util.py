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
Tests for idaes.core.dmf.util module
"""
import datetime
import hashlib
import json
import logging
import os
from pathlib import Path
from tempfile import TemporaryDirectory
import time
from typing import Union

import pytest

#
from idaes.core.dmf import util, resource
from .util import init_logging

__author__ = "Dan Gunter"

init_logging()
_log = logging.getLogger(__name__)
scratch_dir: Union[str, None] = None
scratch_path: Union[Path, None] = None


def setup_module(module):
    global scratch_dir, scratch_path
    scratch_dir = TemporaryDirectory(prefix="idaes.core.dmf_")  # easier to remove later
    scratch_path = Path(scratch_dir.name)


def teardown_module(module):
    pass


@pytest.mark.unit
def test_strlist():
    input = [1, 2, 3]
    output = util.strlist(input, sep="/")
    assert output == "1/2/3"


@pytest.mark.unit
def test_get_file():
    f1 = util.get_file(__file__)
    assert f1 is not None
    f2 = util.get_file(f1)
    assert f2 is not None
    assert f2.name == f1.name


@pytest.mark.unit
def test_import_module():
    m = util.import_module("idaes.core.dmf.util")
    assert m is not None


@pytest.mark.unit
def test_get_module_version():
    m = util.import_module("idaes.core.dmf.util")
    v1 = util.get_module_version(m)
    assert v1 is None
    m.__version__ = "foobar"
    try:
        util.get_module_version(m)
    except ValueError:
        pass
    else:
        assert False, "ValueError expected for {}".format(m.__version__)
    m.__version__ = "1.2.3-alpha"
    util.get_module_version(m)


@pytest.mark.unit
def test_get_module_author():
    m = util.import_module("idaes.core.dmf.util")
    util.get_module_author(m)


@pytest.mark.unit
def test_tempdir():
    with util.TempDir() as newdir:
        pass


@pytest.mark.unit
def test_datetime_timestamp():
    ts = time.time()
    dt = datetime.datetime.fromtimestamp(ts)
    ts1 = util.datetime_timestamp(dt)
    assert pytest.approx(ts, 0.000001) == ts1


@pytest.mark.unit
def test_mkdir_p():
    tmp_dir = scratch_path / "mkdir_p"
    random_str = hashlib.sha1().hexdigest()
    # test absolute
    path = str(tmp_dir / f"{random_str}/idaes.core.dmf/util/mkdir_p")
    util.mkdir_p(path)
    assert os.path.exists(path)
    # test relative
    random_str = hashlib.sha1().hexdigest()
    os.chdir(str(tmp_dir))
    path = f"{random_str}/idaes.core.dmf/util/mkdir_p"
    util.mkdir_p(path)
    assert os.path.exists(path)


@pytest.mark.unit
def test_is_jupyter_notebook():
    nbtext = """
    {
        "cells": [],
        "metadata": {},
        "nbformat": 4,
        "nbformat_minor": 2
    }
    """
    tmp_dir = scratch_path / "is_jupyter_notebook"
    tmp_dir.mkdir()
    f = (tmp_dir / "sample.txt").open("w")
    assert not util.is_jupyter_notebook(f.name)
    f = (tmp_dir / "null.ipynb").open("w")
    assert util.is_jupyter_notebook(f.name, check_contents=False)
    assert not util.is_jupyter_notebook(f.name)
    f = (tmp_dir / "basic.ipynb").open("w")
    f.write(nbtext)
    f.close()
    assert util.is_jupyter_notebook(f.name)
    f = (tmp_dir / "random.ipynb").open("w")
    f.write('{"random": [1,2,3]}')
    f.close()
    assert not util.is_jupyter_notebook(f.name)


@pytest.mark.unit
def test_is_python():
    tmp_dir = scratch_path / "is_python"
    tmp_dir.mkdir()
    f = (tmp_dir / "sample.txt").open("w")
    assert not util.is_python(f.name)
    f = (tmp_dir / "sample.py").open("w")
    assert util.is_python(f.name)


@pytest.mark.unit
def test_is_resource_json():
    tmp_dir = scratch_path / "is_resource_json"
    tmp_dir.mkdir()
    f = (tmp_dir / "sample.txt").open("w")
    assert not util.is_resource_json(f.name)
    f = (tmp_dir / "sample.json").open("w")
    assert not util.is_resource_json(f.name)
    r = resource.Resource()
    json.dump(r.v, f)
    f.close()
    assert util.is_resource_json(f.name)
    assert not util.is_resource_json(f.name, max_bytes=1)
