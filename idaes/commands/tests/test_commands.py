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
def test_examples_cli_noop(runner):
    result = runner.invoke(examples.get_examples, ["--no-install", "--no-download"])
    assert result.exit_code == 0


@pytest.mark.integration()
def test_examples_cli_list(runner):
    result = runner.invoke(examples.get_examples, ["-l"])
    assert result.exit_code == 0


def can_write(path):
    """Test ability to write a file in 'path'.
    Assume the path is a temporary directory, so don't bother cleaning up.
    """
    if not path.exists():
        _log.debug(f"can_write: Creating parent directory '{path}'")
        try:
            path.mkdir()
        except Exception as err:
            _log.warning(f"Failed to create temporary directory '{path}': {err}")
            return False
    test_file = path / "test_file.txt"
    _log.debug(f"can_write: Creating file '{test_file}'")
    try:
        fp = test_file.open("w")
        fp.write("hello")
        fp.close()
    except Exception as err:
        _log.warning(f"Failed to write sample file '{test_file}': {err}")
        return False


@pytest.mark.integration()
def test_examples_cli_download(runner, tempdir):
    # failure with existing dir
    result = runner.invoke(examples.get_examples, ["-d", str(tempdir), "-I"])
    assert result.exit_code == -1


@pytest.mark.integration()
def test_examples_cli_explicit_version(runner, tempdir):
    dirname = str(tempdir / "examples")
    result = runner.invoke(examples.get_examples, ["-d", dirname, "-I", "-V", "1.5.0"])
    assert result.exit_code == 0


@pytest.mark.skip("Erroneously fails when latest ideas version is a pre-release")
@pytest.mark.integration()
def test_examples_cli_default_version(runner, tempdir):
    dirname = str(tempdir / "examples")
    result = runner.invoke(examples.get_examples, ["-d", dirname])
    assert result.exit_code == 0


@pytest.mark.integration()
def test_examples_cli_download_unstable(runner, tempdir):
    if can_write(tempdir):
        dirname = str(tempdir / "examples")
        # unstable version but no --unstable flag
        result = runner.invoke(
            examples.get_examples, ["-d", dirname, "-V", "1.2.3-beta"]
        )
        assert result.exit_code == -1


@pytest.mark.integration()
def test_examples_cli_copy(runner, tempdir):
    dirname = str(tempdir / "examples")
    # local path is bad
    badpath = str(tempdir / "no-such-dir")
    result = runner.invoke(examples.get_examples, ["-d", dirname, "--local", badpath])
    assert result.exit_code == -1
    # local dir exists, no REPO_DIR in it
    src_dir = tempdir / "examples-dev"
    src_dir.mkdir(exist_ok=True)
    result = runner.invoke(
        examples.get_examples, ["-d", dirname, "--local", str(src_dir), "-I"]
    )
    assert result.exit_code == -1
    # local path ok, target is ok
    # src_dir = tempdir / "examples-dev"
    # src_dir.mkdir()
    (src_dir / examples.REPO_DIR).mkdir()
    (src_dir / examples.REPO_DIR / "file.py").open("w")
    result = runner.invoke(
        examples.get_examples, ["-d", dirname, "--local", str(src_dir), "-I"]
    )
    assert result.exit_code == 0
    # make sure the file was copied into the target dir
    assert (Path(dirname) / "file.py").exists()


# non-CLI


@pytest.mark.integration()  # goes out to network
def test_examples_n():
    target_dir = str(uuid.uuid4())  # pick something that won't exist
    retcode = subprocess.call(["idaes", "get-examples", "-N", "-d", target_dir])
    assert retcode != 0  # failure


@pytest.mark.integration()  # goes out to network
def test_examples_list_releases():
    releases = examples.get_releases(True)
    assert len(releases) > 0
    examples.print_releases(releases, True)


@pytest.mark.integration()  # goes out to network
def test_examples_download_bad_version():
    assert pytest.raises(examples.DownloadError, examples.download, Path("."), "1.2.3")


@pytest.mark.integration()
def test_examples_find_python_directories(tmp_path):
    root = tmp_path
    # populate a/c/file.py, a/d/file.py, b/c/file.py, b/d/file.py
    for i in ("a", "b"):
        # create parent
        os.mkdir(str(root / i))
        for j in ("c", "d"):
            pydir = root / i / j
            # create child dir
            os.mkdir(str(pydir))
            # put python file in child
            (pydir / "file.py").open("w")
    # find
    found_dirs = examples.find_python_directories(root)
    # check
    rel_found_dirs = [d.relative_to(root) for d in found_dirs]
    for i in ("a", "b"):
        path_i = Path(i)
        assert path_i in rel_found_dirs
        for j in ("c", "d"):
            assert path_i / j in rel_found_dirs


@pytest.mark.integration()
def test_examples_check_github_response():
    # ok result
    examples.check_github_response([{"result": "ok"}], {})
    # rate limit
    pytest.raises(
        examples.GithubError,
        examples.check_github_response,
        {"message": "API rate limit exceeded,dude"},
        {"X-RateLimit-Reset": "2000000000"},
    )
    # unknown error
    pytest.raises(
        examples.GithubError,
        examples.check_github_response,
        {"message": "Something is seriously wrong"},
        {},
    )
    # unknown, no message - still an error
    pytest.raises(
        examples.GithubError,
        examples.check_github_response,
        {"oscar": "green grouch"},
        {},
    )
    # non-dict: error
    pytest.raises(
        examples.GithubError,
        examples.check_github_response,
        12,
        "hello",
    )


@pytest.mark.nowin32
@pytest.mark.integration
def test_examples_install_src():
    tempdir = create_module_scratch("examples_install_src")
    # monkey patch a random install package name so as not to have
    # any weird side-effects on other tests
    orig_install_pkg = examples.INSTALL_PKG
    examples.INSTALL_PKG = "i" + str(uuid.uuid4()).replace("-", "_")
    _log.debug(f"install_src: curdir={os.curdir}")
    # create fake package
    src_dir = tempdir / "src"
    src_dir.mkdir(exist_ok=True)
    m1_dir = src_dir / "module1"
    m1_dir.mkdir(exist_ok=True)
    (m1_dir / "groot.py").open("w").write("print('I am groot')\n")
    m2_dir = m1_dir / "module1_1"
    m2_dir.mkdir(exist_ok=True)
    (m2_dir / "groot.py").open("w").write("print('I am groot')\n")
    # install it
    examples.install_src("0.0.0", src_dir)
    # patch back the proper install package name
    examples.INSTALL_PKG = orig_install_pkg


@pytest.mark.unit
def test_examples_cleanup(tempdir):
    # put some crap in the temporary dir
    #
    # <tempdir>/
    #   idaes_examples.egg-info/
    #      afile
    #   dist/
    #      idaes_examples-1.0.0-py3.8.egg
    #      idaes_examples-2.0.0-py3.9.egg
    #   <tempdir2>/
    #       afile
    #
    os.chdir(tempdir)
    eggy = Path("idaes_examples.egg-info")
    eggy.mkdir()
    (eggy / "afile").open("w")
    distdir = Path("dist")
    distdir.mkdir()
    (distdir / "idaes_examples-1.0.0-py3.8.egg").open("w")
    (distdir / "idaes_examples-2.0.0-py3.9.egg").open("w")
    name2 = str(uuid.uuid4())
    tempsubdir = Path(name2)
    tempsubdir.mkdir()
    (tempsubdir / "afile").open("w")
    # All of the above should be removed
    # by the cleanup function.
    examples.g_tempdir = tempsubdir
    examples.g_egg = eggy
    examples.clean_up_temporary_files()
    # Check that everything is removed
    # assert not distdir.exists()
    # assert not eggy.exists()
    # assert not tempsubdir.exists()


@pytest.mark.unit
def test_examples_cleanup_nodist(tempdir):
    tempdir = tempdir
    # put some crap in the temporary dir
    #
    # <tempdir>/
    #   idaes_examples.egg-info/
    #      afile
    #   <tempdir2>/
    #       afile
    #
    os.chdir(tempdir)
    eggy = Path("idaes_examples.egg-info")
    eggy.mkdir()
    (eggy / "afile").open("w")
    name2 = str(uuid.uuid4())
    tempsubdir = Path(name2)
    tempsubdir.mkdir()
    (tempsubdir / "afile").open("w")
    # All of the above should be removed
    # by the cleanup function.
    examples.g_tempdir = tempsubdir
    examples.g_egg = eggy
    examples.clean_up_temporary_files()
    # Check that everything is removed
    assert not eggy.exists()
    assert not tempsubdir.exists()


@pytest.mark.unit
def test_examples_cleanup_nodist_noegg(tempdir):
    tempdir = tempdir
    # put some crap in the temporary dir
    #
    # <tempdir>/
    #   <tempdir2>/
    #       afile
    #
    os.chdir(tempdir)
    name2 = str(uuid.uuid4())
    tempsubdir = Path(name2)
    tempsubdir.mkdir()
    (tempsubdir / "afile").open("w")
    # All of the above should be removed
    # by the cleanup function.
    examples.g_tempdir = tempsubdir
    examples.clean_up_temporary_files()
    # Check that everything is removed
    assert not tempsubdir.exists()


@pytest.mark.unit
def test_examples_local(tempdir):
    d = tempdir
    tgt = d / "examples"
    src = d / "src"
    root = d
    # source doesn't exist yet
    pytest.raises(examples.CopyError, examples.copy_contents, tgt, root)
    # source isn't a directory
    src.open("w")
    pytest.raises(examples.CopyError, examples.copy_contents, tgt, root)
    # source exists = OK
    src.unlink()
    src.mkdir()
    (src / "a").open("w")  # create a file
    examples.copy_contents(tgt, root)
    assert (tgt / "a").exists()
    # target now exists, not ok
    pytest.raises(examples.CopyError, examples.copy_contents, tgt, root)
    # done


@pytest.mark.integration()
def test_illegal_dirs(tmp_path):
    root = tmp_path
    # git
    (root / ".git").mkdir()
    try:
        examples.download(root, "1.2.3")
    except examples.DownloadError as err:
        assert "exists" in str(err)
    finally:
        (root / ".git").rmdir()


@pytest.mark.integration()
def test_get_examples_version():
    assert examples.get_examples_version("1.5.0") == "1.5.1"
    assert examples.get_examples_version("foo") == None


test_cell_nb = {
    "cells": [
        {
            "cell_type": "code",
            "execution_count": 2,
            "metadata": {},
            "outputs": [
                {"name": "stdout", "output_type": "stream", "text": ["Hello, world\n"]}
            ],
            "source": ['print("Hello, world")'],
        },
        {
            "cell_type": "code",
            "execution_count": 3,
            "metadata": {"tags": ["test", "remove_cell"]},
            "outputs": [],
            "source": ["assert 2 + 2 == 4"],
        },
    ],
    "metadata": {
        "celltoolbar": "Tags",
        "kernelspec": {
            "display_name": "Python 3",
            "language": "python",
            "name": "python3",
        },
        "language_info": {
            "codemirror_mode": {"name": "ipython", "version": 3},
            "file_extension": ".py",
            "mimetype": "text/x-python",
            "name": "python",
            "nbconvert_exporter": "python",
            "pygments_lexer": "ipython3",
            "version": "3.8.1",
        },
    },
    "nbformat": 4,
    "nbformat_minor": 2,
}


@pytest.fixture
def remove_cells_notebooks(tmp_path):
    # note: first entry must be in root, one entry per subdir, depth search order
    notebooks = [
        tmp_path / d
        for d in (
            "notebook_testing.ipynb",
            Path("a") / "notebook_testing.ipynb",
            Path("a") / "b" / "notebook_testing.ipynb",
        )
    ]
    first = True
    for nb in notebooks:
        if not first:
            nb.parent.mkdir()
        json.dump(test_cell_nb, nb.open("w"))
        first = False
    yield notebooks
    _remove(tmp_path, notebooks)


@pytest.fixture
def remove_cells_notebooks_nosuffix(tmp_path):
    # note: first entry must be in root, one entry per subdir, depth search order
    notebooks = [
        tmp_path / d
        for d in (
            "notebook.ipynb",
            Path("a") / "notebook.ipynb",
            Path("a") / "b" / "notebook.ipynb",
        )
    ]
    first = True
    for nb in notebooks:
        if not first:
            nb.parent.mkdir()
        json.dump(test_cell_nb, nb.open("w"))
        first = False
    yield notebooks
    _remove(tmp_path, notebooks)


def _remove(tmp_path, notebooks):
    for f in tmp_path.rglob("*.ipynb"):
        f.unlink()
    removed = set()
    for nb in reversed(notebooks):
        if nb.parent not in removed:
            nb.parent.rmdir()
            removed.add(nb.parent)


@pytest.mark.unit
def test_find_notebook_files(remove_cells_notebooks):
    root = remove_cells_notebooks[0].parent
    nbfiles = examples.find_notebook_files(root)
    assert len(nbfiles) == len(remove_cells_notebooks)
    for nbfile in nbfiles:
        assert nbfile in remove_cells_notebooks


@pytest.mark.unit
def test_strip_test_cells(remove_cells_notebooks):
    root = remove_cells_notebooks[0].parent
    examples.strip_special_cells(root)
    nbfiles = examples.find_notebook_files(root)
    assert len(nbfiles) == len(remove_cells_notebooks)
    for nbfile in nbfiles:
        with nbfile.open("r") as f:
            nbdata = json.load(f)
            if not nbfile.stem.endswith("_testing"):
                for cell in nbdata["cells"]:
                    cell_meta = cell["metadata"]
                    tags = cell_meta.get("tags", [])
                    assert examples.REMOVE_CELL_TAG not in tags  # tag is gone
            else:
                n = 0
                for cell in nbdata["cells"]:
                    cell_meta = cell["metadata"]
                    tags = cell_meta.get("tags", [])
                    if examples.REMOVE_CELL_TAG in tags:
                        n += 1
                assert n > 0  # tag still there


##################
# get-extensions #
##################


@pytest.mark.unit
def test_get_extensions(runner):
    result = runner.invoke(extensions.get_extensions, ["--no-download"])
    assert result.exit_code == 0


@pytest.mark.unit
def test_print_extensions_version(runner):
    result = runner.invoke(extensions.extensions_version)
    assert result.exit_code == 0


@pytest.mark.unit
def test_get_extensions_plat(runner):
    result = runner.invoke(extensions.bin_platform)
    assert result.exit_code == 0


@pytest.mark.unit
def test_get_extensions_bad_plat(runner):
    result = runner.invoke(extensions.bin_platform, ["--distro", "johns_good_linux42"])
    assert result.exit_code == 0
    assert result.output == "No supported binaries found.\n"


@pytest.mark.unit
def test_extensions_license(runner):
    result = runner.invoke(extensions.extensions_license)
    assert result.exit_code == 0


#################
# convergence  #
################


@pytest.mark.unit
def test_conv_search(runner):
    result = runner.invoke(convergence.convergence_search)
    assert result.exit_code == 0


@pytest.mark.unit
def test_conv_sample(runner):
    fname = os.path.join(idaes.testing_directory, "sample.json")
    result = runner.invoke(
        convergence.convergence_sample,
        [
            "-e",
            "PressureChanger",
            "-N",
            "10",
            "-s",
            fname,
        ],
    )
    assert result.exit_code == 0
    if os.path.exists(fname):
        os.remove(fname)


@pytest.mark.integration
def test_conv_eval(runner):
    fname = os.path.join(idaes.testing_directory, "sample.json")
    fname2 = os.path.join(idaes.testing_directory, "result.json")
    result = runner.invoke(
        convergence.convergence_sample,
        [
            "-e",
            "PressureChanger",
            "-N",
            "10",
            "-s",
            fname,
        ],
    )
    assert result.exit_code == 0
    result = runner.invoke(convergence.convergence_eval, ["-s", fname, "-j", fname2])
    assert result.exit_code == 0
    with open(fname2, "r") as f:
        d = json.load(f)
    assert "inputs" in d
    assert len(d["time_successful"]) == 10
    if os.path.exists(fname):
        os.remove(fname)
    if os.path.exists(fname2):
        os.remove(fname2)


@pytest.mark.unit
def test_conf_display(runner):
    result = runner.invoke(config.config_display)
    assert result.exit_code == 0


@pytest.mark.unit
def test_conf_file_paths(runner):
    result = runner.invoke(config.config_file)
    assert result.exit_code == 0
    result = runner.invoke(config.config_file, ["--global"])
    assert result.exit_code == 0
    result = runner.invoke(config.config_file, ["--local"])
    assert result.exit_code == 0


@pytest.mark.unit
def test_conf_file_paths(runner):
    fname = os.path.join(idaes.testing_directory, "conf_test.json")
    result = runner.invoke(config.config_write, ["--file", fname])
    assert result.exit_code == 0
    with open(fname, "r") as f:
        d = json.load(f)
        assert d["logger_capture_solver"] == True
    if os.path.exists(fname):
        os.remove(fname)


@pytest.mark.unit
def test_conf_set(runner):
    fname = os.path.join(idaes.testing_directory, "conf_test.json")

    def _tst(args):
        result = runner.invoke(
            config.config_set,
            ["logging:loggers:idaes.solver:handlers", "['console']"] + args,
        )
        assert result.exit_code == 0
        with open(fname, "r") as f:
            d = json.load(f)
            assert len(d["logging"]["loggers"]["idaes.solver"]["handlers"]) == 1
            assert d["logging"]["loggers"]["idaes.solver"]["handlers"][0] == "console"
        result = runner.invoke(
            config.config_set,
            ["logging:loggers:idaes.solver:handlers", "'console'", "--del"] + args,
        )
        assert result.exit_code == 0
        result = runner.invoke(
            config.config_set,
            ["logging:loggers:idaes.solver:handlers", "'console'", "--add"] + args,
        )
        assert result.exit_code == 0
        with open(fname, "r") as f:
            d = json.load(f)
            assert len(d["logging"]["loggers"]["idaes.solver"]["handlers"]) == 1
            assert d["logging"]["loggers"]["idaes.solver"]["handlers"][0] == "console"
        result = runner.invoke(
            config.config_set, ["ipopt_l1:options:max_iter", "100"] + args
        )
        assert result.exit_code == 0
        assert "ConfigDict" in str(type(idaes.cfg.ipopt_l1))
        assert "ConfigDict" in str(type(idaes.cfg.ipopt_l1.options))
        assert idaes.cfg.ipopt_l1.options.max_iter == 100
        if os.path.exists(fname):
            os.remove(fname)

    # lbianchi-lbl: these sets of args seem to fail if tensorflow has been imported (#633)
    # TODO these could be refactored as a parametrized pytest.fixture
    if "tensorflow" not in sys.modules:
        _tst(["--global", "--file", fname, "--file_as_global"])
        _tst(["--local", "--file", fname, "--file_as_local"])
    _tst(["--file", fname])


##############
# env info   #
##############


@pytest.mark.unit
def test_env_info1(runner):
    result = runner.invoke(env_info.environment_info)
    assert result.exit_code == 0
