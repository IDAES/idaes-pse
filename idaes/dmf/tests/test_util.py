##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
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
Tests for idaes.dmf.util module
"""
import datetime
import hashlib
import json
import logging
import os
import shutil
import sys
import time

import pytest

#
from idaes.dmf import util, resource
from .util import init_logging, TempDir

if sys.platform.startswith("win"):
    pytest.skip("skipping DMF tests on Windows", allow_module_level=True)

init_logging()
_log = logging.getLogger(__name__)


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
    m = util.import_module("idaes.dmf.util")
    assert m is not None


@pytest.mark.unit
def test_get_module_version():
    m = util.import_module("idaes.dmf.util")
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
    m = util.import_module("idaes.dmf.util")
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
    assert pytest.approx(ts, ts1, 0.000001)


@pytest.mark.unit
def test_mkdir_p():
    random_str = hashlib.sha1().hexdigest()
    # test absolute
    path = "/tmp/{}/idaes/dmf/util/mkdir_p".format(random_str)
    util.mkdir_p(path)
    assert os.path.exists(path)
    shutil.rmtree("/tmp/{}".format(random_str))
    # test relative
    path = "{}/idaes/dmf/util/mkdir_p".format(random_str)
    util.mkdir_p(path)
    assert os.path.exists(path)
    shutil.rmtree(random_str)


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
    with TempDir() as d:
        f = open(os.path.join(d, "sample.txt"), "w")
        assert not util.is_jupyter_notebook(f.name)
        f = open(os.path.join(d, "null.ipynb"), "w")
        assert util.is_jupyter_notebook(f.name, check_contents=False)
        assert not util.is_jupyter_notebook(f.name)
        f = open(os.path.join(d, "basic.ipynb"), "w")
        f.write(nbtext)
        f.close()
        assert util.is_jupyter_notebook(f.name)
        f = open(os.path.join(d, "random.ipynb"), "w")
        f.write('{"random": [1,2,3]}')
        f.close()
        assert not util.is_jupyter_notebook(f.name)


@pytest.mark.unit
def test_is_python():
    with TempDir() as d:
        f = open(os.path.join(d, "sample.txt"), "w")
        assert not util.is_python(f.name)
        f = open(os.path.join(d, "sample.py"), "w")
        assert util.is_python(f.name)


@pytest.mark.unit
def test_is_resource_json():
    with TempDir() as d:
        f = open(os.path.join(d, "sample.txt"), "w")
        assert not util.is_resource_json(f.name)
        f = open(os.path.join(d, "sample.json"), "w")
        assert not util.is_resource_json(f.name)
        r = resource.Resource()
        json.dump(r.v, f)
        f.close()
        assert util.is_resource_json(f.name)
        assert not util.is_resource_json(f.name, max_bytes=1)
