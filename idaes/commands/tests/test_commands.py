#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Tests for idaes.commands
"""
# stdlib
from functools import partial
import json
import logging
import os
from pathlib import Path
from shutil import rmtree
import subprocess
import sys
from typing import Union
import uuid

# third-party
from click.testing import CliRunner
import pytest

# package
from idaes.commands import examples, extensions, convergence, config, env_info, base
from . import create_module_scratch, rmtree_scratch
import idaes

__author__ = "Dan Gunter"


_log = logging.getLogger(__name__)
if os.environ.get("IDAES_TEST_DEBUG", False):
    _log.setLevel(logging.DEBUG)


@pytest.fixture(scope="module")
def runner():
    return CliRunner()


scratch_path: Union[Path, None] = None


def setup_module(module):
    global scratch_path
    scratch_path = create_module_scratch(module.__name__)
    _log.info(f"Using scratch dir: {scratch_path}")


def teardown_module(module):
    _log.info(f"Remove files from: {scratch_path}")
    rmtree_scratch(scratch_path)


@pytest.fixture
def tempdir(request):
    function_name = request.function.__name__[5:]  # remove "test_" prefix
    sub_path = scratch_path / function_name
    if sub_path.exists():
        for f in sub_path.glob("*"):
            rmtree(f)
    else:
        sub_path.mkdir()
    return sub_path


class TestBaseCommand:
    @pytest.fixture
    def run_idaes(self, runner):
        return partial(runner.invoke, base.command_base)

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "flags",
        [
            ["--version"],
            ["--help"],
        ],
        ids=" ".join,
    )
    def test_flags(self, run_idaes, flags):
        result = run_idaes(flags)
        assert result.exit_code == 0


################
# get-examples #
################


@pytest.mark.unit
def test_examples_cli(runner):
    result = runner.invoke(examples.get_examples)
    assert result.exit_code == 0


##################
# get-extensions #
##################
@pytest.mark.integration
def test_get_extensions(runner):
    result = runner.invoke(extensions.get_extensions, ["--no-download"])
    assert result.exit_code == 0


@pytest.mark.integration
def test_print_extensions_version(runner):
    result = runner.invoke(extensions.extensions_version)
    assert result.exit_code == 0


@pytest.mark.integration
def test_get_extensions_plat(runner):
    result = runner.invoke(extensions.bin_platform)
    assert result.exit_code == 0


@pytest.mark.integration
def test_get_extensions_bad_plat(runner):
    result = runner.invoke(extensions.bin_platform, ["--distro", "johns_good_linux42"])
    assert result.exit_code == 0
    assert result.output == "No supported binaries found for johns_good_linux42.\n"


@pytest.mark.integration
def test_extensions_license(runner):
    result = runner.invoke(extensions.extensions_license)
    assert result.exit_code == 0


###########
# config  #
###########


@pytest.mark.unit
def test_conf_display(runner):
    result = runner.invoke(config.config_display)
    assert result.exit_code == 0


@pytest.mark.unit
def test_conf_file_paths(runner):
    result = runner.invoke(config.config_file)
    assert result.exit_code == 0


@pytest.mark.unit
def test_conf_file_paths(runner):
    result = runner.invoke(config.config_write)
    assert result.exit_code == 0


@pytest.mark.skipif(sys.version_info[:2] == (3, 11), reason="Fails on Python 3.11")
@pytest.mark.unit
def test_conf_set(runner):
    result = runner.invoke(config.config_set)
    assert result.exit_code == 0


##############
# env info   #
##############


@pytest.mark.unit
def test_env_info1(runner):
    result = runner.invoke(env_info.environment_info)
    assert result.exit_code == 0
