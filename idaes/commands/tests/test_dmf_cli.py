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
Tests for DMF CLI
"""
# stdlib
import json
import logging
import os
from pathlib import Path
from typing import Union

# third-party
from click.testing import CliRunner
import pytest

# package
from idaes.core.dmf.cli import init, register, info, ls, related, rm, status
from idaes.core.dmf.dmfbase import DMFConfig, DMF
from idaes.core.dmf.workspace import Workspace
from idaes.core.dmf import resource

from . import create_module_scratch, rmtree_scratch

__author__ = "Dan Gunter"

_log = logging.getLogger(__name__)


@pytest.fixture(scope="module")
def runner():
    return CliRunner()


scratch_path: Union[Path, None] = None
dmf_context_num = 1
# on_windows = sys.platform == "win32"  -- maybe needed later


def setup_module(module):
    global scratch_path
    scratch_path = create_module_scratch(module.__name__)


def teardown_module(module):
    rmtree_scratch(scratch_path)


DATAFILE = "foo.txt"


@pytest.fixture()
def dmf_context():
    """Switch DMF context to a random subdir, then switch back when done."""
    global dmf_context_num
    os.chdir(os.path.expanduser("~"))  # make sure we start in HOME
    path = scratch_path / str(dmf_context_num)
    try:
        path.mkdir()
    except:
        pass
    dmf_path = (path / ".dmf").absolute()
    DMFConfig._filename = str(dmf_path)
    origdir = os.getcwd()
    os.chdir(str(path))
    with open(DATAFILE, "w") as fp:
        fp.write("This is some sample data")
    fp.close()
    yield path
    os.unlink(DATAFILE)
    DMFConfig._filename = str(Path("~/.dmf").expanduser())
    os.chdir(origdir)
    dmf_context_num += 1


def create_foo_workspace(runner):
    if (Path("ws") / Workspace.WORKSPACE_CONFIG).exists():
        create_flag = ""
    else:
        create_flag = "--create"
    result = runner.invoke(
        init,
        ["ws", create_flag, "--name", "foo", "--desc", "foo workspace description"],
    )
    assert result.exit_code == 0
    return result


@pytest.mark.unit()
def test_dmf_find(dmf_context, runner):
    filename = DATAFILE
    create_foo_workspace(runner)
    # register some objects in a workspace
    result = runner.invoke(register, [filename])
    id_all = result.output.strip()
    id_4 = id_all[:4]
    #
    result = runner.invoke(info, ["--no-color", id_4], catch_exceptions=False)
    assert result.exit_code == 0
    assert filename in result.output
    #
    result = runner.invoke(
        info, ["--no-color", "--format", "json", id_4], catch_exceptions=False
    )
    assert result.exit_code == 0
    assert filename in result.output
    out = result.output.strip()
    assert out.startswith("{") and out.endswith("}")
    assert '"relations"' in out
    #
    result = runner.invoke(
        info, ["--no-color", "--format", "jsonc", id_4], catch_exceptions=False
    )
    assert result.exit_code == 0
    assert filename in result.output
    out = result.output.strip()
    j = json.loads(out)
    assert len(j["datafiles"]) == 1


@pytest.mark.unit()
def test_dmf_init(dmf_context, runner):
    create_foo_workspace(runner)
    assert (Path("ws") / Workspace.WORKSPACE_CONFIG).exists()
    #
    result = runner.invoke(init, ["ws2", "--create"], input="foo\nfoo desc\n")
    assert result.exit_code == 0
    assert (Path("ws2") / Workspace.WORKSPACE_CONFIG).exists()
    #
    result = runner.invoke(init, ["doesnotexist"])
    assert result.exit_code != 0
    assert "not found" in result.output
    #
    os.mkdir("some_random_directory")
    result = runner.invoke(init, ["some_random_directory"])
    assert result.exit_code != 0
    #
    result = runner.invoke(init, ["ws", "--create"])
    assert result.exit_code != 0
    assert "exists" in result.output


@pytest.mark.unit()
def test_dmf_ls(dmf_context, runner):
    create_foo_workspace(runner)
    # add some files
    files = [f"foo1{n}" for n in range(5)]
    files.append("bar1")
    for f in files:
        with open(f, "w") as fp:
            fp.write("This is some sample data")
        runner.invoke(register, [f])  # add file to DMF
    # regular
    result = runner.invoke(ls, ["--no-color"])
    assert result.exit_code == 0
    output1 = result.output
    result = runner.invoke(ls, ["--no-color"])
    assert result.output == output1
    # sort
    result = runner.invoke(ls, ["--no-color", "-S", "modified"])
    assert result.exit_code == 0
    output1 = result.output
    result = runner.invoke(ls, ["--no-color", "--sort", "modified"])
    assert result.output == output1


# Due to issue with tempfiles
@pytest.mark.component()
def test_dmf_register(dmf_context, runner):
    create_foo_workspace(runner)
    # create CSV file
    filename = "file.csv"
    with open(filename, "w") as fp:
        fp.write("index,time,value\n1,0.1,1.0\n2,0.2,1.3\n")
    # csv
    result = runner.invoke(register, ["file.csv", "--info"], catch_exceptions=False)
    assert result.exit_code == 0
    assert filename in result.output
    assert "version" in result.output
    # csv again / no-unique
    result = runner.invoke(
        register,
        [
            "file.csv",
        ],
        catch_exceptions=False,
    )
    assert result.exit_code != 0
    result = runner.invoke(
        register, ["file.csv", "--no-unique"], catch_exceptions=False
    )
    assert result.exit_code == 0
    # json
    not_json = "notreally.json"
    with open(not_json, "w") as fp:
        fp.write("totally bogus\n")
    result = runner.invoke(register, [not_json], catch_exceptions=False)
    assert result.exit_code == 0
    result = runner.invoke(
        register, [not_json, "--strict", "--no-unique"], catch_exceptions=False
    )
    assert result.exit_code != 0
    # notebook
    not_nb = "my.ipynb"
    with open(not_nb, "w") as fp:
        fp.write("foo\n")
    result = runner.invoke(register, [not_nb, "-t", "notebook"])
    assert result.exit_code != 0  # error: not JSON
    result = runner.invoke(register, [not_nb, "-t", "data"])
    assert result.exit_code == 0
    # relations
    for text_file in "shoebox", "shoes", "closet":
        open(f"{text_file}.txt", "w")
    result = runner.invoke(register, ["shoebox.txt"], catch_exceptions=False)
    assert result.exit_code == 0
    shoebox_id = result.output.strip()
    result = runner.invoke(
        register, ["shoes.txt", "--contained", shoebox_id], catch_exceptions=False
    )
    assert result.exit_code == 0
    shoe_id = result.output.strip()
    result = runner.invoke(info, [shoe_id, "--format", "jsonc"])
    assert result.exit_code == 0
    # relations-2
    for text_file in "shoebox", "shoes", "closet":
        open(f"{text_file}.txt", "w")
    result = runner.invoke(
        register,
        ["closet.txt", "--is-subject", "--contained", shoebox_id],
        catch_exceptions=False,
    )
    assert result.exit_code == 0
    closet_id = result.output.strip()
    result = runner.invoke(info, [shoebox_id, "--format", "jsonc"])
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert len(data["relations"]) == 2
    for rel in data["relations"]:
        if rel["role"] == "subject":
            assert rel["identifier"] == shoe_id
        else:
            assert rel["identifier"] == closet_id


@pytest.mark.unit()
def test_dmf_related(dmf_context, runner):
    create_foo_workspace(runner)
    # add the fully-connected 4 resources
    dmf = DMF()
    rlist = [
        resource.Resource(value={"desc": ltr, "aliases": [ltr], "tags": ["graph"]})
        for ltr in "ABCD"
    ]
    A_id = rlist[0].id  # root resource id, used in testcode
    relation = resource.Predicates.uses
    for r in rlist:
        for r2 in rlist:
            if r is r2:
                continue
            resource.create_relation(r, relation, r2)
    for r in rlist:
        dmf.add(r)
    #
    result = runner.invoke(
        related, [A_id, "--no-unicode", "--no-color"], catch_exceptions=False
    )
    assert result.exit_code == 0
    rlines = result.output.split("\n")
    nrelations = sum(1 for _ in filter(lambda s: resource.Predicates.uses in s, rlines))
    assert nrelations == 12  # 3 blocks of (1 + 3)


@pytest.mark.unit()
def test_dmf_rm(dmf_context, runner):
    create_foo_workspace(runner)
    # add some files named `file{1-5}`
    for i in range(1, 6):
        filename = f"file{i}.txt"
        with open(filename, "w") as fp:
            fp.write(f"file #{i}")
        runner.invoke(register, [filename])
    #
    result = runner.invoke(ls, ["--no-color", "--no-prefix"])
    assert result.exit_code == 0
    output1 = result.output
    output1_lines = output1.split("\n")
    rsrc_id = output1_lines[1].split()[0]  # first token
    result = runner.invoke(rm, [rsrc_id, "--yes", "--no-list"], catch_exceptions=False)
    assert result.exit_code == 0
    result = runner.invoke(ls, ["--no-color", "--no-prefix"])
    assert result.exit_code == 0
    output2 = result.output
    output2_lines = output2.split("\n")
    assert len(output2_lines) == len(output1_lines) - 1
    #
    result = runner.invoke(ls, ["--no-color"])
    assert result.exit_code == 0
    output1 = result.output
    output1_lines = output1.split("\n")
    rsrc_id = output1_lines[1].split()[0]  # first token
    result = runner.invoke(rm, [rsrc_id, "--yes", "--no-list"], catch_exceptions=False)
    assert result.exit_code == 0
    result = runner.invoke(ls, ["--no-color", "--no-prefix"])
    assert result.exit_code == 0
    output2 = result.output
    output2_lines = output2.split("\n")
    assert len(output2_lines) == len(output1_lines) - 1


@pytest.mark.unit()
def test_dmf_status(dmf_context, runner):
    create_foo_workspace(runner)
    #
    result = runner.invoke(status, ["--no-color"])
    assert result.exit_code == 0
    assert "settings" in result.output
    assert "name: foo" in result.output
    #
    result = runner.invoke(status, ["--no-color", "--show", "files"])
    assert result.exit_code == 0
    assert "settings" in result.output
    assert "name: foo" in result.output
    assert "files:" in result.output
    #
    result = runner.invoke(
        status, ["--no-color", "--show", "files", "--show", "htmldocs"]
    )
    assert result.exit_code == 0
    assert "settings" in result.output
    assert "name: foo" in result.output
    assert "html" in result.output
    #
    result = runner.invoke(status, ["--no-color", "-a"])
    assert result.exit_code == 0
    assert "settings" in result.output
    assert "name: foo" in result.output
    assert "html" in result.output
    assert "logging:" in result.output
